"""Cascade-funnel design pipeline for SynBio 2026.

Stages
------
1. Library: combine each parent (canonical 238-aa sfGFP and the previous Top-1
   winner-238) with random subsets of the merged thermostability + brightness
   mutation pools.
2. Hard filters: length / charset / M-start / no stop codons.
3. Exclusion-list filter: O(1) hash lookup against the 135K-entry list.
4. Composite scoring: weighted sum of brightness prior, thermostability prior,
   and (optional) ESM-2 log-likelihood.
5. Diverse Top-6: greedy selection with a Hamming-distance diversity constraint.

Run
---
    python -m src.pipeline --team "YourTeamName" --out designs/submission.csv

The pipeline is deterministic given a `--seed`. CPU-only by default.
"""
from __future__ import annotations
import argparse
import csv
import json
import random
from dataclasses import dataclass, field, asdict
from pathlib import Path

import pandas as pd

from .constants import SFGFP, AVGFP
from .validate import validate_sequence, validate_team_submission
from .exclusion import ExclusionIndex
from .mutate import apply_mutations, random_combinations
from .thermostability import (
    candidate_pool as thermo_pool,
    thermostability_score,
    explain as thermo_explain,
    CANDIDATES as THERMO_CANDIDATES,
)
from .brightness import BrightnessScorer
from .esm_score import get_scorer as get_esm_scorer, available as esm_available
from .ml_brightness import MLBrightnessScorer


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
DATA_GFP = DATA_DIR / "GFP_data.xlsx"
DATA_EXCL = DATA_DIR / "Exclusion_List.csv"


# ---------------------------------------------------------------------------
# Mutation pools
# ---------------------------------------------------------------------------
def brightness_pool(scorer: BrightnessScorer, top_k: int = 2,
                    min_effect: float = 0.03,
                    forbid_aas: tuple = ("C",)) -> list:
    """Curate beneficial single substitutions from the data-driven prior.

    Returns (comp_pos, aa, weight) tuples in competition numbering.
    Excludes Cys mutations because the cell-free E. coli system is reducing
    and disulfide bonds are unreliable for thermostability there.
    """
    options = scorer.per_site_options(top_k=top_k, min_func_rate=0.5, min_support=5)
    pool = []
    for comp_pos, items in options.items():
        for aa, eff in items:
            if aa in forbid_aas:
                continue
            if eff >= min_effect:
                pool.append((comp_pos, aa, float(eff)))
    return pool


def merge_pools(thermo: list, bright: list) -> list:
    """Combine two (pos, aa, weight) pools, deduplicating by (pos, aa)."""
    seen = {}
    for pos, aa, w in thermo + bright:
        key = (pos, aa)
        if key not in seen or w > seen[key][2]:
            seen[key] = (pos, aa, w)
    return list(seen.values())


# ---------------------------------------------------------------------------
# Library generation
# ---------------------------------------------------------------------------
@dataclass
class Candidate:
    seq: str
    parent_name: str
    mutations: list
    template: str = ""           # source template (e.g. "amacGFP") if from candidates CSV
    mutation_str: str = ""       # original mutation annotation (e.g. "I11L;K162A;...")
    n_mutations_orig: int = 0    # mutation count as recorded in the candidates CSV
    brightness: float = 0.0
    thermo: float = 0.0
    esm: float = 0.0
    ml: float = 0.0
    composite: float = 0.0


def load_candidates_csv(path: Path, parents: dict) -> list:
    """Load user-supplied candidate sequences from a CSV.

    The CSV must have a 'Sequence' column. Optional columns: Template, Mutation,
    Mutation_Count, Notes are preserved into the candidate's `mutations` field for
    the design log.

    Each candidate's parent is auto-assigned by sequence length (after the
    2026-05-27 correction, both supported parents are 238 aa; the first 238-aa
    parent in `parents` is used). Sequences whose length doesn't match any
    parent are still loaded but get parent_name='unknown' (scorers that need
    a parent will return 0 for them).
    """
    df = pd.read_csv(path)
    if "Sequence" not in df.columns:
        raise ValueError(f"{path} must have a 'Sequence' column. Got: {list(df.columns)}")
    by_len: dict = {}
    for pname, pseq in parents.items():
        by_len.setdefault(len(pseq), pname)

    out = []
    for i, row in df.iterrows():
        seq = str(row["Sequence"]).strip().upper()
        pname = by_len.get(len(seq), "unknown")
        # Preserve any human-readable mutation annotation if present
        mut_anno = row.get("Mutation") if "Mutation" in df.columns else row.get("Notes", "")
        if pd.isna(mut_anno):
            mut_anno = ""
        mut_str = str(mut_anno)
        mut_list = [t.strip() for t in mut_str.replace(",", ";").split(";") if t.strip()]
        template = str(row.get("Template", "")) if "Template" in df.columns else ""
        if pd.isna(template) or template == "nan":
            template = ""
        try:
            n_mut = int(row.get("Mutation_Count", len(mut_list)))
        except (TypeError, ValueError):
            n_mut = len(mut_list)
        out.append(Candidate(seq=seq, parent_name=pname, mutations=mut_list,
                              template=template, mutation_str=mut_str,
                              n_mutations_orig=n_mut))
    return out


def generate_library(parents: dict, pools: dict, n_per_parent: int = 1500,
                     size_range: tuple = (2, 8), seed: int = 17) -> list:
    """For each parent, generate `n_per_parent` random combinatorial variants."""
    out = []
    rng = random.Random(seed)
    for pname, pseq in parents.items():
        pool = pools[pname]
        seqs = random_combinations(
            pseq, [(p, a) for p, a, _ in pool], n=n_per_parent,
            size_range=size_range, seed=rng.randint(0, 2**31 - 1),
        )
        for s in seqs:
            muts = []
            for i, (a, b) in enumerate(zip(pseq, s), 1):
                if a != b:
                    muts.append(f"{a}{i}{b}")
            out.append(Candidate(seq=s, parent_name=pname, mutations=muts))
    # Always include the unmodified parents themselves
    for pname, pseq in parents.items():
        out.append(Candidate(seq=pseq, parent_name=pname, mutations=[]))
    return out


# ---------------------------------------------------------------------------
# Filtering & scoring
# ---------------------------------------------------------------------------
def hard_filter(cands: list) -> list:
    return [c for c in cands if validate_sequence(c.seq).ok]


def exclusion_filter(cands: list, idx: ExclusionIndex) -> list:
    return [c for c in cands if not idx.contains(c.seq)]


def score_all(cands: list, parents: dict, b_scorer: BrightnessScorer,
              esm_scorer, ml_scorer: MLBrightnessScorer,
              w_b: float, w_t: float, w_e: float, w_ml: float) -> None:
    """Mutates candidates in place, filling all scoring axes + composite.

    Candidates with parent_name='unknown' (no length-matching parent registered)
    get all scores set to 0 — they pass through filters but won't rank highly
    and won't crash the pipeline.
    """
    for c in cands:
        if c.parent_name not in parents:
            c.brightness = 0.0
            c.thermo = 0.0
            c.esm = 0.0
            continue
        parent = parents[c.parent_name]
        c.brightness = b_scorer.score(c.seq, parent=parent)
        c.thermo = thermostability_score(c.seq, parent=parent)
        c.esm = esm_scorer.score(c.seq) if w_e > 0 else 0.0

    # Batched ML scoring (predict all candidates per parent in one call)
    if w_ml > 0 and ml_scorer.available():
        from collections import defaultdict
        by_parent = defaultdict(list)
        for i, c in enumerate(cands):
            if c.parent_name in parents:
                by_parent[c.parent_name].append(i)
            # else: leave c.ml at default 0.0
        for pname, idxs in by_parent.items():
            seqs = [cands[i].seq for i in idxs]
            preds = ml_scorer.score_many(seqs, parent=parents[pname])
            for i, p in zip(idxs, preds):
                cands[i].ml = float(p)

    if not cands:
        return

    def _norm(vals):
        lo, hi = min(vals), max(vals)
        rng = (hi - lo) or 1.0
        return [(v - lo) / rng for v in vals]
    bn = _norm([c.brightness for c in cands])
    tn = _norm([c.thermo for c in cands])
    en = _norm([c.esm for c in cands])
    mn = _norm([c.ml for c in cands])
    for c, b, t, e, m in zip(cands, bn, tn, en, mn):
        c.composite = w_b * b + w_t * t + w_e * e + w_ml * m


# ---------------------------------------------------------------------------
# Diverse Top-6 selection
# ---------------------------------------------------------------------------
def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(1 for x, y in zip(a, b) if x != y)


def diverse_top(cands: list, k: int = 6, min_hamming: int = 4) -> list:
    """Greedy: highest composite first, then enforce Hamming distance to all picks.
    Falls back to relaxing the constraint if we can't fill k slots."""
    cands = sorted(cands, key=lambda c: c.composite, reverse=True)
    picked: list = []
    for c in cands:
        if all(hamming(c.seq, p.seq) >= min_hamming for p in picked):
            picked.append(c)
            if len(picked) >= k:
                return picked
    # Relax constraint if needed
    for c in cands:
        if c not in picked:
            picked.append(c)
            if len(picked) >= k:
                break
    return picked


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def write_submission(picks: list, team: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Team_Name", "Seq_ID", "Sequence"])
        for i, c in enumerate(picks, 1):
            w.writerow([team, i, c.seq])


def write_design_log(picks: list, parents: dict, out_path: Path,
                     meta: dict) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    log = {"meta": meta, "designs": []}
    for i, c in enumerate(picks, 1):
        log["designs"].append({
            "seq_id": i,
            "template": c.template,
            "parent": c.parent_name,
            "length": len(c.seq),
            "n_mutations_reported": c.n_mutations_orig,
            "n_mutations_diff_vs_parent": len(c.mutations),
            "mutation_string": c.mutation_str,
            "mutations": c.mutations,
            "brightness_prior_delta_log10": round(c.brightness, 4),
            "thermostability_prior": round(c.thermo, 4),
            "esm_log_likelihood": round(c.esm, 4),
            "ml_brightness_log10": round(c.ml, 4),
            "composite_score": round(c.composite, 4),
            "thermo_explanation": thermo_explain(c.seq, parent=parents[c.parent_name]) if c.parent_name in parents else [],
            "sequence": c.seq,
        })
    out_path.write_text(json.dumps(log, indent=2))


def write_design_report_csv(picks: list, team: str, out_path: Path) -> None:
    """Write a flat CSV that lists every pick's Template, Mutation, scores —
    useful for the Design Concept Document and team reviews. Distinct from
    submission.csv (which stays in the official 3-column format)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "Seq_ID", "Team_Name", "Template", "Parent",
        "Length", "N_Mutations_reported", "N_Mutations_vs_parent",
        "Mutation", "Bright_Prior", "Thermo_Prior",
        "ESM_PLL", "ML_Pred_log10", "Composite_Score", "Sequence",
    ]
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for i, c in enumerate(picks, 1):
            w.writerow({
                "Seq_ID":                 i,
                "Team_Name":              team,
                "Template":               c.template or "(generated)",
                "Parent":                 c.parent_name,
                "Length":                 len(c.seq),
                "N_Mutations_reported":   c.n_mutations_orig,
                "N_Mutations_vs_parent":  len(c.mutations),
                "Mutation":               c.mutation_str or ";".join(c.mutations),
                "Bright_Prior":           f"{c.brightness:+.4f}",
                "Thermo_Prior":           f"{c.thermo:+.4f}",
                "ESM_PLL":                f"{c.esm:+.4f}",
                "ML_Pred_log10":          f"{c.ml:+.4f}",
                "Composite_Score":        f"{c.composite:.4f}",
                "Sequence":               c.seq,
            })


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--team", default="OurTeam", help="Team name for the submission CSV")
    ap.add_argument("--out", default="designs/submission.csv")
    ap.add_argument("--log", default="designs/design_log.json")
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--n-per-parent", type=int, default=1500)
    ap.add_argument("--w-bright", type=float, default=0.20, help="weight for the prior brightness lookup")
    ap.add_argument("--w-thermo", type=float, default=0.40)
    ap.add_argument("--w-esm",    type=float, default=0.10)
    ap.add_argument("--w-ml",     type=float, default=0.30, help="weight for ML brightness predictor")
    ap.add_argument("--min-hamming", type=int, default=4)
    ap.add_argument("--candidates-csv", default=None,
                    help="Path to a CSV with a 'Sequence' column. When set, "
                         "the pipeline scores YOUR sequences instead of "
                         "generating its own. Combine with --merge to do both.")
    ap.add_argument("--merge", action="store_true",
                    help="Merge user-supplied candidates with generated ones. "
                         "Requires --candidates-csv.")
    args = ap.parse_args()

    print("=" * 64)
    print("SynBio 2026 GFP design pipeline")
    print("=" * 64)

    # 1) Parents: corrected canonical sfGFP (238 aa) and 2025 best (238 aa).
    # Both are 238 aa in canonical avGFP/sfGFP numbering after the 2026-05-27
    # competition correction — no inter-parent position shift required.
    top = pd.read_excel(DATA_GFP, sheet_name="beforetopseqs")
    top2025 = top[top["year"] == 2025]["sequence"].tolist()
    seed_2025 = top2025[0]
    parents = {"sfGFP_238": SFGFP, "winner2025_238": seed_2025}

    # 2) Build mutation pools per parent. Positions are 1-indexed canonical numbering.
    b_scorer = BrightnessScorer()
    bright_pool = brightness_pool(b_scorer, top_k=2, min_effect=0.05)
    thermo_pool_list = thermo_pool()  # canonical 238-numbering

    pools = {
        "sfGFP_238":      merge_pools(thermo_pool_list, bright_pool),
        "winner2025_238": merge_pools(thermo_pool_list, bright_pool),
    }
    print(f"  parents : {list(parents)}")
    print(f"  pool sizes: {{ {', '.join(f'{k}: {len(v)}' for k,v in pools.items())} }}")

    # 3) Build the candidate library
    user_cands = []
    if args.candidates_csv:
        cpath = Path(args.candidates_csv)
        if not cpath.is_absolute():
            cpath = REPO_ROOT / cpath
        user_cands = load_candidates_csv(cpath, parents)
        print(f"  loaded {len(user_cands):,} user candidates from {cpath}")
        # Length distribution + parent assignment
        from collections import Counter
        ln_counts = Counter(len(c.seq) for c in user_cands)
        pn_counts = Counter(c.parent_name for c in user_cands)
        print(f"    length histogram: {dict(ln_counts)}")
        print(f"    parent assignment: {dict(pn_counts)}")

    if args.candidates_csv and not args.merge:
        cands = user_cands
    elif args.candidates_csv and args.merge:
        gen = generate_library(parents, pools, n_per_parent=args.n_per_parent, seed=args.seed)
        cands = user_cands + gen
        print(f"  + {len(gen):,} generated candidates  →  total {len(cands):,}")
    else:
        cands = generate_library(parents, pools, n_per_parent=args.n_per_parent, seed=args.seed)

    print(f"  library size: {len(cands):,}")

    # 4) Hard filters
    cands = hard_filter(cands)
    print(f"  after hard filters: {len(cands):,}")

    # 5) Exclusion-list filter
    print("  loading exclusion index ...")
    excl = ExclusionIndex.from_csv(DATA_EXCL)
    print(f"    excluded sequences: {len(excl):,}")
    cands = exclusion_filter(cands, excl)
    print(f"  after exclusion filter: {len(cands):,}")

    # 6) Score
    ml_scorer = MLBrightnessScorer()
    use_esm = esm_available() and args.w_esm > 0
    use_ml = ml_scorer.available() and args.w_ml > 0
    # Only instantiate the ESM scorer when we actually need it — otherwise
    # importing/loading would trigger a model download (and an SSL error on
    # macs whose Python lacks certifi-installed certificates).
    if use_esm:
        esm_scorer = get_esm_scorer()
    else:
        from .esm_score import _NullScorer
        esm_scorer = _NullScorer()
    if not use_esm:
        args.w_esm = 0.0
    if not use_ml:
        args.w_ml = 0.0
    # Re-normalize remaining weights so they sum to 1
    s = args.w_bright + args.w_thermo + args.w_esm + args.w_ml
    if s > 0:
        args.w_bright /= s
        args.w_thermo /= s
        args.w_esm    /= s
        args.w_ml     /= s
    print(f"  ESM available: {esm_available()}  |  ML available: {ml_scorer.available()}")
    print(f"  weights: B={args.w_bright:.2f}  T={args.w_thermo:.2f}  "
          f"E={args.w_esm:.2f}  ML={args.w_ml:.2f}")
    score_all(cands, parents, b_scorer, esm_scorer, ml_scorer,
              args.w_bright, args.w_thermo, args.w_esm, args.w_ml)

    # 7) Diverse Top-6
    picks = diverse_top(cands, k=6, min_hamming=args.min_hamming)
    print(f"  picked: {len(picks)} diverse top sequences")
    for i, c in enumerate(picks, 1):
        print(f"    [{i}] parent={c.parent_name} mut={len(c.mutations)} "
              f"B={c.brightness:+.3f} T={c.thermo:+.3f} E={c.esm:+.3f} "
              f"ML={c.ml:+.3f} comp={c.composite:.3f} len={len(c.seq)}")

    # 8) Validate the FINAL submission as a whole
    records = [{"Team_Name": args.team, "Seq_ID": i + 1, "Sequence": c.seq}
               for i, c in enumerate(picks)]
    final_check = validate_team_submission(records, args.team)
    if not final_check.ok:
        print("  !! FINAL VALIDATION FAILED:")
        for r in final_check.reasons:
            print("    -", r)
        raise SystemExit(2)
    # Cross-check exclusion list one more time
    for c in picks:
        assert not excl.contains(c.seq), f"selected sequence in exclusion list: {c.seq[:30]}..."

    # 9) Write outputs
    out_path = REPO_ROOT / args.out
    log_path = REPO_ROOT / args.log
    report_path = out_path.with_name("design_report.csv")
    write_submission(picks, args.team, out_path)
    write_design_report_csv(picks, args.team, report_path)
    write_design_log(picks, parents, log_path, meta={
        "seed": args.seed,
        "n_per_parent": args.n_per_parent,
        "weights": {"brightness": args.w_bright, "thermo": args.w_thermo,
                     "esm": args.w_esm, "ml": args.w_ml},
        "esm_available": esm_available(),
        "ml_available":  ml_scorer.available(),
        "ml_model_name": ml_scorer.bundle["model_name"] if ml_scorer.available() else None,
        "library_size_after_filters": len(cands),
    })
    print(f"  → wrote {out_path}")
    print(f"  → wrote {report_path}  (per-pick Template / Mutation / scores)")
    print(f"  → wrote {log_path}")
    print("Done.")


if __name__ == "__main__":
    main()
