"""Score a user-supplied CSV of candidate sequences against the trained model.

Workflow
--------
You provide a CSV (default name `test.csv`) with at minimum a `Sequence` column.
Optionally include `Seq_ID`, `Notes`, and `Brightness` (the measured log10
brightness if you have it; used for evaluation metrics only).

The script:
  1. Validates each sequence (length 220-250, M-start, AA20 only, no stop).
  2. Checks each against the 135 414-entry exclusion list.
  3. Computes mutations vs the canonical 238-aa sfGFP (post-2026-05-27
     competition correction). No position shift is needed.
  4. Predicts brightness with the trained ML model.
  5. Computes the thermostability prior from literature mutations.
  6. Computes the data-driven brightness lookup prior.
  7. Combines them into a composite score (same weights as the design pipeline).
  8. If a Brightness column was present, reports per-sequence error and
     overall RMSE / R².
  9. Writes a ranked output CSV.

Usage
-----
From the project root:

    python -m src.predict --input test.csv --out designs/test_predictions.csv

Optional flags:

    --w-bright F   weight for brightness lookup prior (default 0.20)
    --w-thermo F   weight for thermostability prior   (default 0.40)
    --w-ml     F   weight for ML brightness predictor (default 0.40)
    --parent {sfGFP,avGFP,winner2025}
                   parent for mutation diff (default sfGFP)
"""
from __future__ import annotations
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import SFGFP, AVGFP, AA_SET, MIN_LEN, MAX_LEN
from .validate import validate_sequence
from .exclusion import ExclusionIndex
from .thermostability import (
    thermostability_score, explain as thermo_explain, CANDIDATES as THERMO,
)
from .brightness import BrightnessScorer
from .ml_brightness import MLBrightnessScorer

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
DATA_EXCL = DATA_DIR / "Exclusion_List.csv"
PREV_TOP_PATH = DATA_DIR / "GFP_data.xlsx"


def load_parent(name: str) -> tuple:
    if name == "sfGFP":
        return SFGFP, "sfGFP_238"
    if name == "avGFP":
        return AVGFP, "avGFP_238"
    if name == "winner2025":
        # First 2025 entry from beforetopseqs sheet
        df = pd.read_excel(PREV_TOP_PATH, sheet_name="beforetopseqs")
        seq = df[df["year"] == 2025]["sequence"].iloc[0]
        return seq, "winner2025_238"
    raise ValueError(f"Unknown parent: {name}")


def diff_mutations(seq: str, parent: str) -> list:
    """Return list like ['F46L','I167T','V206K'] in canonical 238-aa numbering."""
    if len(seq) != len(parent):
        return [f"<length-mismatch:{len(seq)}vs{len(parent)}>"]
    out = []
    for i, (a, b) in enumerate(zip(parent, seq), 1):
        if a != b:
            out.append(f"{a}{i}{b}")
    return out


def normalize(v: np.ndarray) -> np.ndarray:
    lo, hi = float(v.min()), float(v.max())
    rng = (hi - lo) or 1.0
    return (v - lo) / rng


def main():
    ap = argparse.ArgumentParser(description="Score a CSV of candidate sequences.")
    ap.add_argument("--input", default="test.csv",
                    help="Path to test CSV with at least a `Sequence` column.")
    ap.add_argument("--out", default="designs/test_predictions.csv")
    ap.add_argument("--parent", choices=["sfGFP", "avGFP", "winner2025"], default="sfGFP")
    ap.add_argument("--w-bright", type=float, default=0.20)
    ap.add_argument("--w-thermo", type=float, default=0.40)
    ap.add_argument("--w-ml",     type=float, default=0.40)
    ap.add_argument("--allow-excluded", action="store_true",
                    help="Score even sequences in the exclusion list. Use this for "
                         "model EVALUATION; never for submission candidates.")
    args = ap.parse_args()

    in_path = Path(args.input)
    if not in_path.is_absolute():
        in_path = ROOT / in_path
    if not in_path.exists():
        print(f"[ERROR] Input file not found: {in_path}")
        print("       Pass --input PATH/TO/test.csv (path relative to project root).")
        sys.exit(2)

    print(f"=" * 70)
    print(f"  Scoring candidates from {in_path}")
    print(f"=" * 70)

    df = pd.read_csv(in_path)
    if "Sequence" not in df.columns:
        print("[ERROR] CSV must have a 'Sequence' column.")
        print(f"       Found columns: {list(df.columns)}")
        sys.exit(2)
    print(f"  loaded {len(df)} rows | columns: {list(df.columns)}")

    parent_seq, parent_name = load_parent(args.parent)
    print(f"  parent: {parent_name} ({len(parent_seq)} aa)")

    # ---- 1) Validate + check exclusion list -------------------------------
    print("  loading exclusion index ...")
    excl = ExclusionIndex.from_csv(DATA_EXCL)
    print(f"    {len(excl):,} excluded sequences")

    n_valid, n_excluded = 0, 0
    rows = []
    for i, r in df.iterrows():
        seq = str(r["Sequence"]).strip().upper()
        v = validate_sequence(seq)
        in_excl = excl.contains(seq) if v.ok else False
        muts = diff_mutations(seq, parent_seq) if v.ok else []
        rows.append({
            "Seq_ID":      r.get("Seq_ID", i + 1),
            "Sequence":    seq,
            "Length":      len(seq),
            "Valid":       v.ok,
            "Validation":  "; ".join(v.reasons) if v.reasons else "OK",
            "In_Exclusion": in_excl,
            "N_Mutations": len(muts),
            "Mutations":   ",".join(muts),
        })
        score_eligible = v.ok and (args.allow_excluded or not in_excl)
        if score_eligible:
            n_valid += 1
        if in_excl:
            n_excluded += 1
    if args.allow_excluded:
        print(f"  valid: {n_valid}/{len(df)}  |  in exclusion list: {n_excluded} "
              f"(scoring anyway because --allow-excluded was set)")
    else:
        print(f"  valid+novel: {n_valid}/{len(df)}  |  in exclusion list: {n_excluded}")

    # ---- 2) Score the eligible ones ---------------------------------------
    valid_idx = [i for i, r in enumerate(rows)
                 if r["Valid"] and (args.allow_excluded or not r["In_Exclusion"])]
    valid_seqs = [rows[i]["Sequence"] for i in valid_idx]
    if len(valid_seqs) == 0:
        print(f"\n[ERROR] No sequences passed validation. Nothing to score.")
        print(f"        If you meant to evaluate against the training data, re-run with --allow-excluded.")
        sys.exit(2)

    b_scorer = BrightnessScorer()
    ml_scorer = MLBrightnessScorer()

    bright_arr = np.array([b_scorer.score(s, parent=parent_seq) for s in valid_seqs],
                          dtype=np.float32)
    thermo_arr = np.array([thermostability_score(s, parent=parent_seq) for s in valid_seqs],
                          dtype=np.float32)

    if ml_scorer.available():
        ml_arr = ml_scorer.score_many(valid_seqs, parent=parent_seq)
        ml_model_name = ml_scorer.bundle["model_name"]
    else:
        ml_arr = np.zeros(len(valid_seqs), dtype=np.float32)
        ml_model_name = None
        print("  [WARN] ML model not found at src/_brightness_model.pkl. "
              "Run `python -c 'from src.ml_brightness import train_and_select; train_and_select()'` "
              "or `python run.py --retrain` first.")

    print(f"  ML model: {ml_model_name}")
    print(f"  brightness lookup: range [{bright_arr.min():+.3f}, {bright_arr.max():+.3f}]")
    print(f"  thermo prior     : range [{thermo_arr.min():+.3f}, {thermo_arr.max():+.3f}]")
    print(f"  ML prediction    : range [{ml_arr.min():+.3f}, {ml_arr.max():+.3f}]")

    # Re-normalize weights to account for missing components
    w_b, w_t, w_ml = args.w_bright, args.w_thermo, args.w_ml
    if not ml_scorer.available():
        w_ml = 0.0
    s = w_b + w_t + w_ml
    if s > 0:
        w_b, w_t, w_ml = w_b / s, w_t / s, w_ml / s
    print(f"  weights: B={w_b:.2f}  T={w_t:.2f}  ML={w_ml:.2f}")

    # Normalize each axis (over the valid set) and combine
    if len(valid_idx) > 0:
        bn = normalize(bright_arr)
        tn = normalize(thermo_arr)
        mn = normalize(ml_arr) if ml_scorer.available() else np.zeros_like(ml_arr)
        composite = w_b * bn + w_t * tn + w_ml * mn
    else:
        composite = np.array([])

    # Splat scores back into the rows (invalid rows get 0 / NaN)
    score_lookup = {
        valid_idx[i]: {
            "Bright_Prior":  float(bright_arr[i]),
            "Thermo_Prior":  float(thermo_arr[i]),
            "ML_Pred_log10": float(ml_arr[i]),
            "Composite":     float(composite[i]),
        }
        for i in range(len(valid_idx))
    }
    for idx, r in enumerate(rows):
        s = score_lookup.get(idx, {})
        r["Bright_Prior"]  = s.get("Bright_Prior", float("nan"))
        r["Thermo_Prior"]  = s.get("Thermo_Prior", float("nan"))
        r["ML_Pred_log10"] = s.get("ML_Pred_log10", float("nan"))
        r["Composite"]     = s.get("Composite", float("nan"))

    # ---- 3) (Optional) evaluate against ground truth ---------------------
    have_truth = "Brightness" in df.columns
    rmse = r2 = None
    if have_truth and len(valid_idx) > 0:
        truth = pd.to_numeric(df["Brightness"], errors="coerce").to_numpy()
        for idx, r in enumerate(rows):
            r["Brightness_truth"] = float(truth[idx]) if not np.isnan(truth[idx]) else float("nan")
        # Compute metrics over valid+truth-present subset
        mask = np.array([
            (i in score_lookup) and (not np.isnan(rows[i].get("Brightness_truth", float("nan"))))
            for i in range(len(rows))
        ])
        if mask.sum() >= 2:
            y_true = np.array([rows[i]["Brightness_truth"] for i in range(len(rows)) if mask[i]])
            y_pred = np.array([rows[i]["ML_Pred_log10"]    for i in range(len(rows)) if mask[i]])
            rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
            ss_res = float(np.sum((y_true - y_pred) ** 2))
            ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    # ---- 4) Write output ranked by composite -----------------------------
    out_rows = sorted(rows, key=lambda r: (-(r.get("Composite") or float("-inf")),
                                            r.get("Seq_ID", 0)))
    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = ROOT / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ["Seq_ID", "Length", "Valid", "Validation", "In_Exclusion",
                  "N_Mutations", "Mutations",
                  "Bright_Prior", "Thermo_Prior", "ML_Pred_log10", "Composite",
                  "Sequence"]
    if have_truth:
        fieldnames.insert(-1, "Brightness_truth")

    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    # ---- 5) Console summary ----------------------------------------------
    print()
    print("Top 10 by composite score:")
    print(f"  {'rank':>4}  {'comp':>6}  {'ML':>6}  {'thermo':>7}  {'mut':>4}  {'valid':>5}  Seq_ID")
    for k, r in enumerate(out_rows[:10], 1):
        c = r.get("Composite")
        c = f"{c:+.3f}" if c is not None and not np.isnan(c) else "  n/a"
        ml = r.get("ML_Pred_log10")
        ml = f"{ml:+.3f}" if ml is not None and not np.isnan(ml) else "  n/a"
        th = r.get("Thermo_Prior")
        th = f"{th:+.3f}" if th is not None and not np.isnan(th) else "   n/a"
        flag = "OK" if r["Valid"] and not r["In_Exclusion"] else "FAIL"
        print(f"  {k:>4}  {c:>6}  {ml:>6}  {th:>7}  {r['N_Mutations']:>4}  {flag:>5}  {r['Seq_ID']}")

    if rmse is not None:
        print()
        print(f"Evaluation against ground-truth Brightness column:")
        print(f"  n_evaluated = {int(mask.sum())}")
        print(f"  RMSE = {rmse:.3f}")
        print(f"  R²   = {r2:+.3f}")

    print()
    print(f"Wrote → {out_path}")


if __name__ == "__main__":
    main()
