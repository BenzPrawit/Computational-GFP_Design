#!/usr/bin/env python3
"""One-command terminal runner for the SynBio 2026 GFP design pipeline.

Usage from a Mac/Linux terminal:

    cd /Users/prawitthitayanuwat/Desktop/My-project/SynBio2026
    python3 -m venv .venv && source .venv/bin/activate     # one-time setup
    pip install -r requirements.txt scikit-learn xgboost   # one-time setup
    python run.py --team "YourTeamName"                    # the daily command

What it does:
  1. Trains the ML brightness predictor (RF/GBR/MLP/XGBoost with K-fold CV)
     unless a cached model is already present and you didn't pass --retrain.
  2. Runs the cascade-funnel design pipeline (literature thermo + data-driven
     brightness + ML rerank + exclusion-list filter + diverse Top-6).
  3. Validates the final submission against every competition rule.
  4. Prints a summary and the path to designs/submission.csv.

Common flags:

    --team NAME            team name written into submission.csv  (REQUIRED before submitting)
    --mode {quick,full}    quick=5K subsample, full=141K          (default: quick)
    --retrain              force retraining even if model is cached
    --no-train             skip training, use existing src/_brightness_model.pkl
    --seed N               reproducibility seed                    (default: 17)
    --gridsearch           run hyperparameter grid search          (slow)
    --w-thermo F           weight for thermostability prior        (default: 0.40)
    --w-ml F               weight for ML brightness prediction     (default: 0.40)
    --w-bright F           weight for brightness lookup prior      (default: 0.20)
    --n-per-parent N       library size per parent                 (default: 2500)

Examples:

    # First time, quick test
    python run.py --team "MyTeam"

    # Final submission run with full data and grid search
    python run.py --team "MyTeam" --mode full --gridsearch --retrain
"""
from __future__ import annotations
import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MODEL_PATH = ROOT / "src" / "_brightness_model.pkl"


def banner(text: str) -> None:
    print()
    print("=" * 70)
    print(f"  {text}")
    print("=" * 70)


def run(cmd: list, check: bool = True) -> int:
    print(f"$ {' '.join(cmd)}")
    return subprocess.run(cmd, check=check).returncode


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    ap.add_argument("--team", required=False, default="OurTeam",
                    help='Team name written into submission.csv. REQUIRED before final submission.')
    ap.add_argument("--mode", choices=["quick", "full"], default="quick",
                    help='quick=fast (5K rows), full=best quality (141K rows). Default: quick.')
    ap.add_argument("--retrain", action="store_true",
                    help="Force retraining even if a cached model exists.")
    ap.add_argument("--no-train", action="store_true",
                    help="Skip training entirely (requires existing model).")
    ap.add_argument("--gridsearch", action="store_true",
                    help="Run a small hyperparameter grid search before final fit.")
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--n-per-parent", type=int, default=2500)
    ap.add_argument("--w-bright", type=float, default=0.20)
    ap.add_argument("--w-thermo", type=float, default=0.40)
    ap.add_argument("--w-ml",     type=float, default=0.40)
    ap.add_argument("--w-esm",    type=float, default=0.0)
    ap.add_argument("--out",      default="designs/submission.csv")
    ap.add_argument("--log",      default="designs/design_log.json")
    ap.add_argument("--candidates-csv", default=None,
                    help="Path to a CSV with a 'Sequence' column. The pipeline "
                         "will score and rank YOUR sequences. Add --merge to also "
                         "include pipeline-generated sequences.")
    ap.add_argument("--merge", action="store_true",
                    help="Combine user candidates with pipeline-generated ones.")
    args = ap.parse_args()

    # Helpful warning for the most common forgotten step
    if args.team == "OurTeam":
        print("[WARN] You're using the placeholder team name 'OurTeam'.")
        print("       Pass --team \"YourActualTeamName\" for the real submission.\n")

    # ---- 1) Train (or skip) ------------------------------------------------
    need_train = (not args.no_train) and (args.retrain or not MODEL_PATH.exists())
    if need_train:
        banner("Step 1/2 — train ML brightness model (dual: all_data + bright_focused)")
        if args.mode == "quick":
            subsample, k = 5_000, 5
        else:
            subsample, k = None, 10

        train_code = f"""
from src.ml_brightness import train_dual_model
train_dual_model(family='avGFP', subsample={subsample!r}, k={k})
"""
        t0 = time.time()
        rc = subprocess.run(
            [sys.executable, "-c", train_code], cwd=ROOT,
        ).returncode
        if rc != 0:
            print("[ERROR] training failed.")
            sys.exit(rc)
        print(f"\n[OK] training finished in {time.time()-t0:.1f}s — model saved to {MODEL_PATH}")
    else:
        if MODEL_PATH.exists():
            print(f"[skip] training — using cached model at {MODEL_PATH}")
            print("       (pass --retrain to rebuild it)")
        else:
            print("[ERROR] --no-train was set but no cached model exists at",
                  MODEL_PATH)
            sys.exit(2)

    # ---- 2) Design ---------------------------------------------------------
    banner("Step 2/2 — run cascade-funnel design pipeline")
    cmd = [
        sys.executable, "-m", "src.pipeline",
        "--team", args.team,
        "--seed", str(args.seed),
        "--n-per-parent", str(args.n_per_parent),
        "--w-bright", str(args.w_bright),
        "--w-thermo", str(args.w_thermo),
        "--w-ml",     str(args.w_ml),
        "--w-esm",    str(args.w_esm),
        "--out", args.out,
        "--log", args.log,
    ]
    if args.candidates_csv:
        cmd += ["--candidates-csv", args.candidates_csv]
    if args.merge:
        cmd += ["--merge"]
    rc = run(cmd, check=False)
    if rc != 0:
        print("[ERROR] design pipeline failed.")
        sys.exit(rc)

    # ---- 3) Summary --------------------------------------------------------
    banner("Done")
    sub_path = ROOT / args.out
    log_path = ROOT / args.log
    print(f"  submission   → {sub_path}")
    print(f"  design log   → {log_path}")
    print()
    print("Quick view of the submission CSV:")
    if sub_path.exists():
        for ln in sub_path.read_text().splitlines()[:8]:
            # Truncate long sequence rows for terminal readability
            cells = ln.split(",", 2)
            if len(cells) == 3 and len(cells[2]) > 60:
                cells[2] = cells[2][:60] + f"... ({len(cells[2])} aa)"
            print("  " + ",".join(cells))
    print()
    print("Next steps:")
    print(f"  • Open {sub_path} to inspect the 6 final sequences")
    print(f"  • Open {log_path} for per-sequence provenance (mutations, scores)")
    print('  • Re-run with --mode full --gridsearch for the final submission')


if __name__ == "__main__":
    main()
