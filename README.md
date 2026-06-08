# SynBio 2026 — Computational GFP Design (ML → AlphaFold → AMBER MD)

A reproducible computational pipeline for designing and prioritising green
fluorescent protein (GFP) variants for the **SynBio 2026 Challenge**. The
project combines a machine-learning sequence–function model with a structural
and physics-based validation cascade (AlphaFold 3 → AMBER molecular
dynamics → high-temperature thermostability analysis) to recommend a short list
of candidate sequences.

> **Status — computational predictions only.** Every candidate, score, and
> ranking in this repository is computational. None has been experimentally
> validated. Wet-lab expression and characterisation are required before any
> claim about brightness or thermostability can be made.

---

## Project overview

The objective of the SynBio 2026 Challenge is to design GFP variants that
maximise **initial brightness** while retaining **thermostability** under a
defined heat-stress assay, subject to strict length, format, and
exclusion-list constraints. Rather than de novo generation, this project starts
from established avGFP/sfGFP-like backbones and an internal mutational library,
then filters candidates through a cascade funnel and validates the survivors
structurally and dynamically.

### Workflow

```
Stage 1  ML model training & validation        (src/, notebooks/)
Stage 2  Selection of candidate GFP sequences   (results/tables/)
Stage 3  AlphaFold 3 structure prediction       (structures/alphafold3/)
Stage 4  AMBER MD system setup                   (md_simulation/, structures/prepared/)
Stage 5  High-temperature MD stress test         (md_simulation/)
Stage 6  Trajectory analysis & ranking           (md_simulation/MD-Analysis/, results/)
```

### Research objectives

1. Train and honestly validate an ML model for GFP sequence → brightness, with explicit checks against data leakage.
2. Generate and rank a candidate library, accounting for model uncertainty and applicability domain.
3. Assess structural plausibility of the top candidates with AlphaFold 3.
4. Compare relative structural stability of candidates vs. a reference GFP using AMBER MD under a high-temperature stress protocol.
5. Produce a transparent, reproducible record suitable for a report, thesis, or manuscript.

---

## Directory structure

```
github-upload/
├── README.md                 ← this file
├── LICENSE                    ← MIT
├── CITATION.cff               ← how to cite this repository
├── requirements.txt           ← pip dependencies
├── environment.yml            ← conda environment
├── .gitignore / .gitattributes
│
├── data/
│   ├── raw/                   ← official competition inputs (GFP_data.xlsx, AA seqs, brief, template)
│   ├── processed/             ← derived inputs (test_candidates.csv)
│   └── external/              ← pointers to large/restricted data NOT stored here
│
├── notebooks/                 ← exploratory & end-to-end ML notebooks
│
├── src/                       ← ML pipeline package (validate, mutate, brightness, thermostability, pipeline, …)
│
├── models/                    ← trained model artifact (best_model.pkl)
│
├── structures/
│   ├── reference/             ← reference crystal (2B3P),
│   ├── alphafold3/            ← AF3 predicted models (.cif) + confidence summaries (gfp_1…gfp_6)
│   └── prepared/              ← per-candidate prepared PDBs for MD (GFP-0…GFP-6)
│
├── md_simulation/             ← AMBER minimisation/heating/equilibration/production input decks
│   └── MD-Analysis/           ← CPPTRAJ trajectory-analysis scripts (RMSD, RMSF, Rg, H-bonds, contacts)
│
├── docs/                      ← AMBER workflow guides, AlphaFold terms, project skill spec, REFERENCES.md
│
└── tests/                     ← lightweight reproducibility checks
```

A short note on deviations from the generic data-science template: this project
contains **two coupled computational pipelines** (an ML pipeline that produces
candidate sequences, and a structural/MD pipeline that validates them), so two
extra top-level folders are used. `structures/` holds 3D models that are neither
raw competition data nor final results, and `md_simulation/` holds the AMBER
input decks that drive the simulations.

---

## Installation

Python 3.10+ is recommended.

**With pip:**

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

**With conda:**

```bash
conda env create -f environment.yml
conda activate synbio2026-gfp
```

The structural and MD stages additionally require external scientific software
that is **not** installed by pip/conda: AlphaFold 3 (or AlphaFold Server output)
and AmberTools/AMBER. See `docs/` for the guides.

---

## Usage

### ML pipeline (Stages 1–2)

The pipeline package lives in `src/`. From the repository root:

```bash
# Train the brightness model and generate ranked candidates
python -m src.pipeline --team "YourTeam" --seed 17 \
    --out results/tables/submission.csv \
    --log results/tables/design_log.json
```

`src/run.py` is a one-command convenience runner that trains the model (if no
cached model is present) and runs the design pipeline. Run `python src/run.py --help`
for all flags (`--mode {quick,full}`, `--retrain`, `--seed`, scoring weights, …).

> The original `run.py` / `run.sh` helper scripts assumed a `designs/` output
> folder and a hard-coded macOS path. In this repository, outputs are written to
> `results/tables/`. Update output paths accordingly when you run them.

### Structural & MD validation (Stages 3–6)

These stages are documented step-by-step in:

- `docs/GFP_AMBER_MD_Guideline.md` — project-specific AMBER MD recipe
- `docs/AMBER-Guide.md` — general AMBER reference

The AMBER input decks in `md_simulation/` follow the order:
`min* → heat → eq1…eq4 → md` (production), with the high-temperature stress run
targeting **352.15 K (79 °C)**. Trajectory analysis scripts (RMSD, RMSF, Rg,
H-bonds, contacts) are in `md_simulation/MD-Analysis/`.

---

## Example workflow

```bash
# 1. Set up environment
conda env create -f environment.yml && conda activate synbio2026-gfp

# 2. Train + design (Stages 1–2)
python -m src.pipeline --team "YourTeam" --out results/tables/submission.csv

# 3. Predict AF3 structures for the 6 candidates (external; see docs)
#    → place models under structures/alphafold3/gfp_N/

# 4. Prepare systems and run AMBER MD at 352.15 K (external; see docs)
#    → use md_simulation/ decks; outputs analysed with md_simulation/MD-Analysis/

# 5. Rank candidates on multiple stability metrics (not RMSD alone)
```

## Expected outputs

| Stage | Key output | Location |
|---|---|---|
| ML model | CV report, held-out report | `results/tables/ml_cv_report.json`, `strict_family_holdout_report.json` |
| Candidates | Final 6 sequences (submission format) | `results/tables/submission.csv`, `designed_sequences.fasta` |
| Provenance | Per-sequence design log | `results/tables/design_log.json` |
| AlphaFold | Predicted models + confidence | `structures/alphafold3/gfp_*/` |
| MD | Trajectories & stability metrics | generated locally (not stored — see `.gitignore`) |
| Reporting | Rationale, validation plan, decks | `results/reports/` |

The current `results/tables/submission.csv` contains 6 candidate sequences
(team `NHLG-2`), each a 238-aa avGFP/sfGFP-like backbone.

---

## Reproducibility notes

- Set the `--seed` flag for deterministic ML runs.
- The MD stage depends on external software versions; record AmberTools,
  force field, and water model in your methods (templates in `docs/`).
- Large generated artifacts (MD trajectories, `.dat` analysis files, AF3 PAE
  matrices and MSAs) are intentionally excluded — see `.gitignore` and
  `data/external/README.md`.
- One target-temperature inconsistency exists in the source material: the
  competition assay is described at **72 °C**, while the MD stress protocol uses
  **79 °C / 352.15 K**. These are different conditions (experimental assay vs.
  in-silico stress test); confirm and state both explicitly in any write-up.

---

## Citation

If you use this repository, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff). Published references that informed the design
(superfolder GFP, TGP, StayGold, mBaoJin, the avGFP local fitness landscape) are
listed in [`docs/REFERENCES.md`](docs/REFERENCES.md). Those papers are
copyrighted and are **not** redistributed here.

## License

Released under the [MIT License](LICENSE).

## Contact

Prawit Thitayanuwat, Mahidol University
ORCID: [0009-0009-7209-541X](https://orcid.org/0009-0009-7209-541X)
Email: prawit.tht@student.mahidol.ac.th (institutional) · pwttynwt.8@gmail.com (personal)
