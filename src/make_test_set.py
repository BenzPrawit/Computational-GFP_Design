"""Build a held-out test set with measured ground-truth brightness.

Sources two kinds of test rows so you can measure both ranking and absolute
accuracy of the predictor:

  1. **Random hold-out from GFP_data.xlsx** (avGFP family)
     A deterministic 10% slice of the rows, converted from mutation-string
     format ('A109D:N145D:...') to full-length sequences vs the avGFP parent.
     IMPORTANT: re-train the model with these rows excluded if you want a true
     held-out evaluation (use --exclude-from-training output).

  2. **Previous Top sequences** (the `beforetopseqs` sheet)
     Known-good designs from 2024 + 2025. They have no measured brightness in
     the dataset, but they are *implicitly high-brightness* since they won
     previous rounds. Useful as positive controls.

Usage
-----
    python -m src.make_test_set --out test_eval.csv --frac 0.10 --seed 42

The output CSV has Seq_ID, Sequence, Brightness, Notes — exactly the shape
that `src/predict.py` consumes.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import pandas as pd

from .constants import AVGFP

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
DATA_PATH = DATA_DIR / "GFP_data.xlsx"


def parse_mut_string(s: str) -> list:
    if not isinstance(s, str) or s.upper() == "WT" or not s.strip():
        return []
    out = []
    for tok in s.split(":"):
        tok = tok.strip()
        if len(tok) < 3: continue
        try:
            out.append((tok[0], int(tok[1:-1]), tok[-1]))
        except (ValueError, IndexError):
            pass
    return out


def apply_mutations_avgfp(parent: str, muts: list) -> str:
    s = list(parent)
    for ref, pos, alt in muts:
        i = pos - 1
        if 0 <= i < len(s):
            s[i] = alt
    return "".join(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="test_eval.csv")
    ap.add_argument("--frac", type=float, default=0.10,
                    help="Fraction of avGFP rows to set aside as test set.")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--max-test-rows", type=int, default=2000,
                    help="Cap test set size to keep predict.py fast.")
    ap.add_argument("--include-prev-top", action="store_true",
                    help="Append the 20 previous Top-10 sequences as positive controls.")
    ap.add_argument("--exclude-from-training", default=None,
                    help="If set, write a CSV of avGFP indices to EXCLUDE from training.")
    args = ap.parse_args()

    print(f"Loading {DATA_PATH} ...")
    df = pd.read_excel(DATA_PATH, sheet_name="brightness")
    av = df[df["GFP type"] == "avGFP"].copy().reset_index(drop=True)
    av = av[av["aaMutations"].astype(str).str.upper() != "WT"]
    print(f"  avGFP rows: {len(av):,}")

    # Deterministic 10% slice via seeded sampling
    test = av.sample(frac=args.frac, random_state=args.seed).head(args.max_test_rows)
    train = av.drop(test.index)
    print(f"  train: {len(train):,}  test: {len(test):,}")

    rows = []
    skipped = 0
    for sid, (_, r) in enumerate(test.iterrows(), 1):
        muts = parse_mut_string(str(r["aaMutations"]))
        if not muts:
            skipped += 1
            continue
        seq = apply_mutations_avgfp(AVGFP, muts)
        # Skip if any mutation position is out of range (rare data noise)
        if any(p < 1 or p > len(AVGFP) for _, p, _ in muts):
            skipped += 1
            continue
        rows.append({
            "Seq_ID":     sid,
            "Sequence":   seq,
            "Brightness": float(r["Brightness"]),
            "Notes":      f"avGFP+{r['aaMutations']}",
        })
    print(f"  built {len(rows)} test rows ({skipped} skipped)")

    if args.include_prev_top:
        top = pd.read_excel(DATA_PATH, sheet_name="beforetopseqs")
        for _, r in top.iterrows():
            rows.append({
                "Seq_ID":     len(rows) + 1,
                "Sequence":   r["sequence"],
                "Brightness": "",
                "Notes":      f"prev_top_{r['year']}",
            })
        print(f"  + appended {len(top)} previous Top sequences (no brightness label)")

    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = ROOT / out_path
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["Seq_ID", "Sequence", "Brightness", "Notes"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote → {out_path}  ({len(rows)} rows)")

    if args.exclude_from_training:
        excl_path = Path(args.exclude_from_training)
        if not excl_path.is_absolute():
            excl_path = ROOT / excl_path
        # Save the raw mutation strings of the held-out rows so the trainer can
        # filter them out explicitly during data loading.
        excl_path.write_text("\n".join(test["aaMutations"].astype(str).tolist()))
        print(f"Wrote held-out mutation strings → {excl_path}")
        print(f"(Train with these rows removed for a true held-out evaluation.)")


if __name__ == "__main__":
    main()
