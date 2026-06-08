# External / large / restricted data — pointers

Some inputs are **not stored in this repository** because they are large,
copyrighted, or regenerable. This file records what they are and where to get
them, so the project remains reproducible.

| Item | Why excluded | How to obtain / regenerate | Original location |
|---|---|---|---|
| `Exclusion_List.csv` (~32 MB) | Large; competition-provided | Download from the SynBio 2026 competition portal, then place in `data/external/`. Track with Git LFS if you must commit it. | `ml_part/data/Exclusion_List.csv` |
| Published reference papers (PDFs) | Copyrighted — must not be redistributed | See full citations and DOIs in [`../../docs/REFERENCES.md`](../../docs/REFERENCES.md) | `ml_part/referencepaper/` |
| `Amber-2026.pdf` (~20 MB AMBER manual) | Large; redistributable terms restricted | Download from https://ambermd.org/ | `qm_md_part/Amber-2026.pdf` |
| AMBER MD trajectories & `.dat` analysis files (~18 MB+) | Regenerable simulation output | Regenerate by running the AMBER pipeline in `config/amber/` | `qm_md_part/Protein/GFP-0/Out/` |
| AlphaFold 3 `*_full_data_*.json` (PAE) + `msas/`, `templates/` | Large; regenerable from the job request | Re-run AlphaFold Server/AF3 with the kept `*_job_request.json` | `qm_md_part/Protein/AlphaFold3/gfp_*/` |

## Notes

- The pipeline's exclusion-identity check requires `Exclusion_List.csv`. Place it
  in `data/external/` before running the full design pipeline.
- AlphaFold output is subject to its own terms — see
  [`../../docs/AlphaFold3_terms_of_use.md`](../../docs/AlphaFold3_terms_of_use.md).
- All excluded paths are listed in the repository `.gitignore`.
