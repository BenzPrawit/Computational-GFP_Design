"""ML brightness predictor — adapted from the user's QSAR workflow.

Pipeline (mirrors the QSAR_ML.ipynb structure):
  1. Featurize sequences (sequence-based descriptors instead of molecular ones)
  2. Train RF / GradientBoosting / XGBoost / MLP regressors
  3. 10-fold cross-validation to pick the best model by RMSE / R²
  4. Save the winning model + featurizer for the cascade-funnel pipeline to load
  5. Optional SHAP analysis for interpretability

Key adaptation note
-------------------
Your QSAR work was small-data (~28 compounds × many molecular descriptors), so
VIF for collinearity made sense. Here we have ~141K mutations × engineered
sequence descriptors — the *opposite* regime. We replace VIF with feature
selection by frequency (only one-hot mutations seen ≥ N times) and rely on
tree-model regularization to handle correlation.
"""
from __future__ import annotations
import json
import time
import pickle
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
DATA_PATH = DATA_DIR / "GFP_data.xlsx"
MODEL_PATH = REPO_ROOT / "src" / "_brightness_model.pkl"
CV_REPORT_PATH = REPO_ROOT / "reports" / "ml_cv_report.json"

# Physicochemical AA properties (Kyte-Doolittle hydrophobicity, volume,
# charge at pH 7). Used to compute "delta" features per mutation.
AA_HYDRO = {  # Kyte-Doolittle
    "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
    "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
}
AA_VOL = {  # side-chain volume Å^3 (Zamyatnin)
    "A":  88, "R": 173, "N": 114, "D": 111, "C": 108,
    "Q": 143, "E": 138, "G":  60, "H": 153, "I": 166,
    "L": 166, "K": 168, "M": 162, "F": 189, "P": 112,
    "S":  89, "T": 116, "W": 227, "Y": 193, "V": 140,
}
AA_CHARGE = {a: 0 for a in "ACDEFGHIKLMNPQRSTVWY"}
AA_CHARGE.update({"D": -1, "E": -1, "K": +1, "R": +1, "H": +0.1})

WT_LABEL = "WT"


def _parse_mut_string(s: str) -> list:
    if not isinstance(s, str) or s.upper() == WT_LABEL or not s.strip():
        return []
    out = []
    for tok in s.split(":"):
        tok = tok.strip()
        if len(tok) < 3:
            continue
        try:
            out.append((tok[0], int(tok[1:-1]), tok[-1]))
        except (ValueError, IndexError):
            pass
    return out


# ---------------------------------------------------------------------------
# Featurizer
# ---------------------------------------------------------------------------
FAMILIES = ("avGFP", "amacGFP", "cgreGFP", "ppluGFP")


class MutationFeaturizer:
    """Convert a mutation set into a fixed-length feature vector.

    Features:
      [0]   number of mutations
      [1]   sum of |Δhydrophobicity|
      [2]   sum of  Δhydrophobicity   (signed)
      [3]   sum of |Δvolume|
      [4]   sum of  Δvolume           (signed)
      [5]   sum of  Δcharge
      [6]   number of mutations to charged residues (D,E,K,R)
      [7]   number of mutations to hydrophobic residues (A,V,L,I,M,F,W)
      [8]   number of mutations to Pro
      [9]   number of mutations to Gly
      [10..13]  family one-hot (avGFP / amacGFP / cgreGFP / ppluGFP)
                — only used when use_family=True; zero otherwise
      [14..N]  one-hot of frequent mutations (pos|alt seen ≥ min_count times,
               keyed by FAMILY when use_family=True so cross-family rows don't
               collide on the same position)
    """

    BASE_GLOBAL_DIM = 10
    FAMILY_DIM = len(FAMILIES)

    def __init__(self, min_count: int = 5, max_features: int = 2000,
                 use_family: bool = False):
        self.min_count = min_count
        self.max_features = max_features
        self.use_family = use_family
        self.mut_index: dict = {}
        self.feature_names: list = []

    @property
    def GLOBAL_DIM(self) -> int:
        return self.BASE_GLOBAL_DIM + (self.FAMILY_DIM if self.use_family else 0)

    @staticmethod
    def _mut_key(family: str | None, pos: int, alt: str) -> tuple:
        return (family, pos, alt)

    def fit(self, mutation_sets: list, families: list | None = None) -> "MutationFeaturizer":
        """fit on either raw mutation lists OR (family-aware) (family, muts) pairs."""
        counter = Counter()
        for i, muts in enumerate(mutation_sets):
            fam = families[i] if (self.use_family and families is not None) else None
            for ref, pos, alt in muts:
                counter[self._mut_key(fam, pos, alt)] += 1
        kept = [k for k, c in counter.most_common(self.max_features) if c >= self.min_count]
        self.mut_index = {k: i for i, k in enumerate(kept)}
        names = ["n_mut", "abs_dHydro", "dHydro", "abs_dVol", "dVol",
                 "dCharge", "n_to_charged", "n_to_hydroph", "n_to_Pro", "n_to_Gly"]
        if self.use_family:
            names += [f"family_{f}" for f in FAMILIES]
        for k in kept:
            fam, p, a = k
            tag = f"{fam}_" if fam else ""
            names.append(f"mut_{tag}{p}{a}")
        self.feature_names = names
        return self

    def transform(self, mutation_sets: list, families: list | None = None) -> np.ndarray:
        n = len(mutation_sets)
        d = self.GLOBAL_DIM + len(self.mut_index)
        X = np.zeros((n, d), dtype=np.float32)
        fam_idx = {f: i for i, f in enumerate(FAMILIES)}
        for i, muts in enumerate(mutation_sets):
            X[i, 0] = len(muts)
            for ref, pos, alt in muts:
                if ref in AA_HYDRO and alt in AA_HYDRO:
                    dh = AA_HYDRO[alt] - AA_HYDRO[ref]
                    dv = AA_VOL[alt] - AA_VOL[ref]
                    dq = AA_CHARGE[alt] - AA_CHARGE[ref]
                    X[i, 1] += abs(dh); X[i, 2] += dh
                    X[i, 3] += abs(dv); X[i, 4] += dv
                    X[i, 5] += dq
                    if alt in "DEKR":      X[i, 6] += 1
                    if alt in "AVLIMFW":   X[i, 7] += 1
                    if alt == "P":         X[i, 8] += 1
                    if alt == "G":         X[i, 9] += 1
                fam = families[i] if (self.use_family and families is not None) else None
                idx = self.mut_index.get(self._mut_key(fam, pos, alt))
                if idx is not None:
                    X[i, self.GLOBAL_DIM + idx] = 1.0
            if self.use_family and families is not None:
                fi = fam_idx.get(families[i])
                if fi is not None:
                    X[i, self.BASE_GLOBAL_DIM + fi] = 1.0
        return X

    def fit_transform(self, mutation_sets: list, families: list | None = None) -> np.ndarray:
        return self.fit(mutation_sets, families).transform(mutation_sets, families)


# ---------------------------------------------------------------------------
# Data loader
# ---------------------------------------------------------------------------
def load_dataset(family: str = "avGFP", subsample: int | None = None,
                 random_state: int = 0,
                 exclude_file: str | Path | None = None) -> tuple:
    """Return (mutation_sets, y, wt_brightness, families_list) for the chosen family.

    Pass `family='all'` to use ALL four families (avGFP + amacGFP + cgreGFP +
    ppluGFP). `wt_brightness` is then the avGFP WT (used as the reference
    baseline because sfGFP is closest to avGFP), and `families_list` carries
    the per-row family label so the featurizer can encode it.

    `exclude_file`: optional path to mutation strings to drop from training.
    """
    df = pd.read_excel(DATA_PATH, sheet_name="brightness")
    if family != "all":
        df = df[df["GFP type"] == family].copy()
    else:
        df = df.copy()
    wt_b = float(df[(df["GFP type"] == "avGFP")
                     & (df["aaMutations"].astype(str).str.upper() == WT_LABEL)
                    ]["Brightness"].mean())
    df = df[df["aaMutations"].astype(str).str.upper() != WT_LABEL]

    if exclude_file is not None:
        ep = Path(exclude_file)
        if not ep.is_absolute():
            ep = REPO_ROOT / ep
        if not ep.exists():
            print(f"[WARN] exclude_file {ep} not found; nothing excluded.")
        else:
            excluded = set(line.strip() for line in ep.read_text().splitlines() if line.strip())
            n_before = len(df)
            df = df[~df["aaMutations"].astype(str).isin(excluded)]
            print(f"  excluded {n_before - len(df):,} rows from training "
                  f"(per {ep.name}, {len(excluded):,} mutation strings listed)")

    if subsample is not None and subsample < len(df):
        df = df.sample(n=subsample, random_state=random_state)
    mutation_sets = [_parse_mut_string(s) for s in df["aaMutations"].astype(str)]
    y = df["Brightness"].to_numpy(dtype=np.float32)
    families_list = df["GFP type"].astype(str).tolist()
    return mutation_sets, y, wt_b, families_list


# ---------------------------------------------------------------------------
# Models — same zoo as the QSAR notebook
# ---------------------------------------------------------------------------
def get_models(seed: int = 7) -> dict:
    """Return a dict of (name -> sklearn-compatible regressor) for the model zoo.

    Includes the QSAR-classic four (RF / GBR / MLP / XGB) plus two modern additions
    ported from the ChatGPT dual-model notebook:
      - ExtraTreesRegressor   — randomized RF cousin; often a small accuracy bump
      - HistGradientBoostingRegressor — sklearn's modern fast booster (LightGBM-style)

    MLP is wrapped in StandardScaler to avoid scale-sensitivity issues.
    """
    from sklearn.ensemble import (
        RandomForestRegressor, ExtraTreesRegressor,
        GradientBoostingRegressor, HistGradientBoostingRegressor,
    )
    from sklearn.neural_network import MLPRegressor
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    try:
        from xgboost import XGBRegressor
        have_xgb = True
    except ImportError:
        have_xgb = False

    models = {
        "RandomForest": RandomForestRegressor(
            n_estimators=220, max_depth=24, min_samples_leaf=3,
            max_features=0.68, n_jobs=-1, random_state=seed,
        ),
        "ExtraTrees": ExtraTreesRegressor(
            n_estimators=260, max_features=0.72, min_samples_leaf=2,
            n_jobs=-1, random_state=seed,
        ),
        "GradientBoosting": GradientBoostingRegressor(
            n_estimators=260, max_depth=5, learning_rate=0.055,
            subsample=0.85, random_state=seed,
        ),
        "HistGradientBoosting": HistGradientBoostingRegressor(
            max_iter=350, learning_rate=0.055, l2_regularization=0.02,
            random_state=seed, max_leaf_nodes=31,
        ),
        "MLP": Pipeline([
            ("scale", StandardScaler(with_mean=False)),  # with_mean=False handles sparse safely
            ("mlp",   MLPRegressor(
                hidden_layer_sizes=(96, 48), activation="relu",
                solver="adam", alpha=1e-4, max_iter=120,
                early_stopping=True, random_state=seed,
            )),
        ]),
    }
    if have_xgb:
        models["XGBoost"] = XGBRegressor(
            n_estimators=400, max_depth=7, learning_rate=0.055,
            subsample=0.85, colsample_bytree=0.85, reg_lambda=1.2,
            objective="reg:squarederror", tree_method="hist",
            n_jobs=-1, random_state=seed, verbosity=0,
        )
    return models


# ---------------------------------------------------------------------------
# K-fold cross-validation (mirrors notebook Part 2 / 3 / 4 / 5)
# ---------------------------------------------------------------------------
def kfold_cv(model, X: np.ndarray, y: np.ndarray, k: int = 10,
             random_state: int = 5) -> dict:
    from sklearn.model_selection import KFold
    from sklearn.metrics import mean_squared_error, r2_score
    kf = KFold(n_splits=k, shuffle=True, random_state=random_state)
    rmse_list, r2_list = [], []
    for fold, (tr, te) in enumerate(kf.split(X), 1):
        model.fit(X[tr], y[tr])
        yp = model.predict(X[te])
        rmse = float(np.sqrt(mean_squared_error(y[te], yp)))
        r2 = float(r2_score(y[te], yp))
        rmse_list.append(rmse)
        r2_list.append(r2)
    return {
        "rmse_mean": float(np.mean(rmse_list)),
        "rmse_std":  float(np.std(rmse_list)),
        "r2_mean":   float(np.mean(r2_list)),
        "r2_std":    float(np.std(r2_list)),
        "rmse_per_fold": rmse_list,
        "r2_per_fold":   r2_list,
    }


# ---------------------------------------------------------------------------
# Train + select + save (mirrors notebook "fix parameters" + final fit cells)
# ---------------------------------------------------------------------------
def train_dual_model(family: str = "avGFP", subsample: int | None = None,
                     bright_threshold: float = 2.5,
                     min_count: int = 5, max_features: int = 2000,
                     k: int = 10, seed: int = 7,
                     exclude_file: str | Path | None = None) -> dict:
    """Train TWO models on the same featurization, addressing the brightness floor:

      * `all_data_model`     — trained on every row (good for "is this dead?" signal)
      * `bright_focused_model` — trained only on functional rows
        (brightness >= `bright_threshold`, default 2.5)

    Both are saved into the same bundle. `MLBrightnessScorer.score_many()` averages
    their predictions when both are present, weighting bright-focused 60% / all-data 40%.

    Use this in place of `train_and_select()` for the production model.
    """
    print(f"[1/5] Loading {family} data ...")
    mut_sets, y, wt_b, fams = load_dataset(family=family, subsample=subsample,
                                            random_state=seed, exclude_file=exclude_file)
    use_family = (family == "all")
    print(f"      n samples = {len(y):,}   WT log10 brightness = {wt_b:.3f}")
    if use_family:
        from collections import Counter as _C
        print(f"      family counts: {dict(_C(fams))}")

    print(f"[2/5] Featurizing (min_count={min_count}, max_features={max_features}, "
          f"use_family={use_family}) ...")
    feat = MutationFeaturizer(min_count=min_count, max_features=max_features,
                               use_family=use_family)
    X = feat.fit_transform(mut_sets, families=fams if use_family else None)
    print(f"      X shape = {X.shape}")

    print(f"[3/5] {k}-fold CV across {len(get_models(seed=seed))}-model zoo on ALL data ...")
    models = get_models(seed=seed)
    cv_results = {}
    timings = {}
    for name, model in models.items():
        t0 = time.time()
        from sklearn.base import clone
        cv = kfold_cv(clone(model), X, y, k=k, random_state=5)
        cv["time_seconds"] = time.time() - t0
        cv_results[name] = cv
        print(f"      {name:>22s}  RMSE={cv['rmse_mean']:.3f}±{cv['rmse_std']:.3f}  "
              f"R²={cv['r2_mean']:+.3f}  ({cv['time_seconds']:.1f}s)")

    best_name = min(cv_results, key=lambda n: cv_results[n]["rmse_mean"])
    print(f"      best model: {best_name}")

    print(f"[4/5] Training all_data_model + bright_focused_model with {best_name} ...")
    bright_mask = y >= bright_threshold
    n_bright = int(bright_mask.sum())
    print(f"      bright rows (brightness >= {bright_threshold}): {n_bright:,}/{len(y):,} "
          f"({100*n_bright/len(y):.1f}%)")

    from sklearn.base import clone
    all_model = clone(get_models(seed=seed)[best_name])
    all_model.fit(X, y)

    bright_model = clone(get_models(seed=seed)[best_name])
    bright_model.fit(X[bright_mask], y[bright_mask])

    # CV the bright-focused model independently
    cv_bright = kfold_cv(clone(get_models(seed=seed)[best_name]),
                          X[bright_mask], y[bright_mask], k=min(k, max(2, n_bright // 500)),
                          random_state=5)
    print(f"      bright-focused CV: RMSE={cv_bright['rmse_mean']:.3f}  R²={cv_bright['r2_mean']:+.3f}")

    print("[5/5] Saving bundle ...")
    bundle = {
        "model":             all_model,           # legacy alias for all_data_model
        "all_data_model":    all_model,
        "bright_focused_model": bright_model,
        "bright_threshold":  bright_threshold,
        "featurizer":        feat,
        "model_name":        best_name,
        "wt_brightness":     wt_b,
        "family":            family,
        "report": {
            "family": family, "n_samples": int(len(y)), "wt_brightness": wt_b,
            "feature_dim": int(X.shape[1]), "n_global": feat.GLOBAL_DIM,
            "n_mutation_onehots": len(feat.mut_index), "k": k,
            "best_model": best_name,
            "models": cv_results,
            "bright_focused_cv": cv_bright,
            "bright_threshold": bright_threshold,
            "n_bright": n_bright,
            "exclude_file": str(exclude_file) if exclude_file else None,
        },
    }
    MODEL_PATH.parent.mkdir(parents=True, exist_ok=True)
    with MODEL_PATH.open("wb") as f:
        pickle.dump(bundle, f)
    CV_REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    CV_REPORT_PATH.write_text(json.dumps(bundle["report"], indent=2))
    print(f"      saved → {MODEL_PATH}")
    print(f"      cv report → {CV_REPORT_PATH}")
    return bundle


def train_and_select(family: str = "avGFP", subsample: int | None = 30_000,
                     min_count: int = 5, max_features: int = 2000,
                     k: int = 10, seed: int = 7) -> dict:
    print(f"[1/4] Loading {family} data ...")
    mut_sets, y, wt_b, _fams = load_dataset(family=family, subsample=subsample,
                                             random_state=seed)
    print(f"      n samples = {len(y):,}   WT log10 brightness = {wt_b:.3f}")

    print(f"[2/4] Fitting featurizer (min_count={min_count}, max_features={max_features}) ...")
    feat = MutationFeaturizer(min_count=min_count, max_features=max_features)
    X = feat.fit_transform(mut_sets)
    print(f"      X shape = {X.shape}   |   "
          f"global features = {feat.GLOBAL_DIM}, "
          f"one-hot mutation features = {len(feat.mut_index)}")

    print(f"[3/4] {k}-fold cross-validation across model zoo ...")
    models = get_models(seed=seed)
    report = {"family": family, "n_samples": int(len(y)), "wt_brightness": wt_b,
              "feature_dim": int(X.shape[1]), "n_global": feat.GLOBAL_DIM,
              "n_mutation_onehots": len(feat.mut_index), "k": k, "models": {}}
    timings = {}
    for name, model in models.items():
        t0 = time.time()
        cv = kfold_cv(model, X, y, k=k, random_state=5)
        elapsed = time.time() - t0
        timings[name] = elapsed
        report["models"][name] = cv
        print(f"      {name:>17s}  RMSE={cv['rmse_mean']:.3f}±{cv['rmse_std']:.3f}  "
              f"R²={cv['r2_mean']:+.3f}±{cv['r2_std']:.3f}  ({elapsed:.1f}s)")

    # Pick the model with the lowest mean RMSE
    best_name = min(report["models"], key=lambda n: report["models"][n]["rmse_mean"])
    print(f"[4/4] Best model = {best_name}. Refitting on full dataset ...")
    best = get_models(seed=seed)[best_name]
    best.fit(X, y)

    bundle = {
        "model": best,
        "featurizer": feat,
        "model_name": best_name,
        "wt_brightness": wt_b,
        "family": family,
        "report": report,
    }
    MODEL_PATH.parent.mkdir(parents=True, exist_ok=True)
    with MODEL_PATH.open("wb") as f:
        pickle.dump(bundle, f)
    CV_REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    CV_REPORT_PATH.write_text(json.dumps(report, indent=2))
    print(f"      saved → {MODEL_PATH}")
    print(f"      cv report → {CV_REPORT_PATH}")
    return bundle


# ---------------------------------------------------------------------------
# Inference helper used by the design pipeline
# ---------------------------------------------------------------------------
class MLBrightnessScorer:
    """Loads the trained bundle and predicts log10 brightness for candidate sequences.

    Positions are 1-indexed canonical avGFP/sfGFP numbering (1..238). After the
    2026-05-27 competition correction the sfGFP reference is canonical 238 aa,
    so no position shift is needed. The legacy 239-aa branch is retained for
    backwards compatibility if a caller passes the OLD erroneous parent.
    """

    def __init__(self, bundle_path: Path = MODEL_PATH):
        if not bundle_path.exists():
            self.bundle = None
            return
        with bundle_path.open("rb") as f:
            self.bundle = pickle.load(f)

    def available(self) -> bool:
        return self.bundle is not None

    def _seq_to_mut_set(self, seq: str, parent: str) -> list:
        if len(seq) != len(parent):
            return []
        out = []
        is_239_legacy = (len(parent) == 239)  # only true with OLD erroneous parent
        for pos, (a, b) in enumerate(zip(parent, seq), 1):
            if a == b:
                continue
            if is_239_legacy:
                if pos <= 171:    canon = pos
                elif pos == 172:  continue
                else:             canon = pos - 1
            else:
                canon = pos  # canonical 238-aa parent → no shift
            out.append((a, canon, b))
        return out

    def score_many(self, seqs: list, parent: str,
                   bright_weight: float = 0.6) -> np.ndarray:
        """Predict brightness for many sequences.

        If the bundle contains both `all_data_model` and `bright_focused_model`
        (dual-model bundle from `train_dual_model`), the prediction is a weighted
        average: `bright_weight * bright_pred + (1-bright_weight) * all_pred`.
        Otherwise falls back to whatever single model is in the bundle.
        """
        if not self.available():
            return np.zeros(len(seqs), dtype=np.float32)
        mut_sets = [self._seq_to_mut_set(s, parent) for s in seqs]
        X = self.bundle["featurizer"].transform(mut_sets)

        bright = self.bundle.get("bright_focused_model")
        all_m  = self.bundle.get("all_data_model") or self.bundle.get("model")
        if bright is not None and all_m is not None:
            return (bright_weight * bright.predict(X)
                    + (1 - bright_weight) * all_m.predict(X)).astype(np.float32)
        # Single-model legacy bundle
        return self.bundle.get("model", all_m).predict(X).astype(np.float32)


# ---------------------------------------------------------------------------
# Optional SHAP analysis (mirrors notebook Part 3 SHAP section)
# ---------------------------------------------------------------------------
def shap_top_features(bundle_path: Path = MODEL_PATH, n_background: int = 100,
                       n_explain: int = 200, top: int = 25) -> list:
    try:
        import shap
    except ImportError:
        print("(install shap to enable: pip install shap)")
        return []
    with bundle_path.open("rb") as f:
        bundle = pickle.load(f)
    feat = bundle["featurizer"]
    model = bundle["model"]
    mut_sets, y, _, _fams = load_dataset(family=bundle["family"], subsample=n_explain + n_background)
    X = feat.transform(mut_sets)
    bg, X_exp = X[:n_background], X[n_background:n_background + n_explain]
    explainer = shap.TreeExplainer(model, data=bg)
    sv = explainer.shap_values(X_exp)
    mean_abs = np.abs(sv).mean(axis=0)
    order = np.argsort(mean_abs)[::-1][:top]
    return [(feat.feature_names[i], float(mean_abs[i])) for i in order]


if __name__ == "__main__":
    train_and_select(family="avGFP", subsample=30_000)
