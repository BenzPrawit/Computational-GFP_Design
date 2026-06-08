# GFP AMBER MD Simulation Guideline

**Project:** SynBio 2026 — GFP All-Around Champion
**Stage:** Post-ML candidate validation (AlphaFold → **AMBER MD** → thermostability)
**Folder:** `md_simulation/` (input decks) and `md_simulation/MD-Analysis/` (CPPTRAJ)
**Plan note (2026-06-07):** The QM stage has been removed. Chromophore parameters now come from a **published Amber-compatible parameter set** (or AM1-BCC charges via `antechamber`, which needs no QM engine). Path references to `qm_md_part/`, `Command/`, and `Command/CPPTRAJ/` below correspond to `md_simulation/` and `md_simulation/MD-Analysis/` in this repository.
**Author:** Prawit Thitayanuwat
**Date:** 2026-05-24

> **How this guide was built**
> This guide is built **directly on the files already present in `qm_md_part/`** — `AMBER-Guide.md`, `Amber-2026.pdf`, and the `Command/` and `Command/CPPTRAJ/` scripts (`addH.in`, `addWater.in`, `min_H.in`, `minWAT.in`, `minALL.in`, `min1.i`–`min4.i`, `heat.i`, `eq1.i`, `eq2.i`, `md.i`, `hbond.in`, `mmpbsa.in`, `*.job`, `min_heat.qsub`, `runMD.qsub`, and the CPPTRAJ scripts).
>
> **Important context:** Those files were originally written for a **protein–small-molecule-ligand complex** (302 residues; ligand at residue 302; manual disulfide bonds; `ZN` in restraint masks). GFP is a **single-chain monomer (~239 residues for competition sfGFP) with a covalently-attached chromophore (CRO at positions 65–67 after cyclization) and no zinc or disulfides**. Every section below flags the **GFP-specific adaptation** needed.
>
> **Scientific caution (applies everywhere):** AMBER MD can support interpretation of structural stability under the simulated conditions. It does **not** prove experimental thermostability or fluorescence. Use language such as *"suggests relative stability under the simulation conditions"* — never *"is thermostable"*.

---

## Table of contents

1. [Project Overview](#1-project-overview)
2. [Required Input Files](#2-required-input-files)
3. [GFP Structure Preparation](#3-gfp-structure-preparation)
4. [AMBER System Building (`tleap`)](#4-amber-system-building-tleap)
5. [Minimization Workflow](#5-minimization-workflow)
6. [Heating and Equilibration](#6-heating-and-equilibration)
7. [Production MD Simulation](#7-production-md-simulation)
8. [Thermostability Testing (79 °C / 352.15 K)](#8-thermostability-testing-79c--35215-k)
9. [CPPTRAJ Analysis](#9-cpptraj-analysis)
10. [Quality Control and Troubleshooting](#10-quality-control-and-troubleshooting)
11. [Final Workflow Summary](#11-final-workflow-summary)

---

## 1. Project Overview

### Purpose of the GFP MD simulation

To assess, in silico, whether each of the 6 ML-selected GFP variants in `ml_part/outputs_notebook/submission.csv` retains a plausible folded structure at ambient temperature (300 K) and under a heat-stress condition (target **79 °C = 352.15 K** per the user-specified protocol; note the competition assay temperature in `ml_part/README.md` is 72 °C — pick one as primary and document it).

### How this connects to the larger GFP design workflow

```
ML candidate selection (ml_part/)
        ↓ submission.csv (6 sequences)
AlphaFold 3 structure prediction
        ↓ <candidate>_AF3.pdb
Chromophore parameterization (published Amber-compatible set, or AM1-BCC via antechamber — no QM)
        ↓ chromophore.mol2 + chromophore.frcmod
AMBER MD (this guide, md_simulation/)
        ↓ trajectories at 300 K and 352.15 K
CPPTRAJ analysis & cross-candidate ranking
        ↓
Final report (with explicit "needs experimental validation" caveat)
```

### Why AMBER

- AMBER's pmemd.cuda gives fast, well-validated classical MD with modern protein force fields (ff19SB) and explicit-solvent water models (OPC). The existing scripts in `Command/` are written for AMBER and use this combination.
- AMBER supports clean integration of nonstandard residues (the cyclic GFP chromophore) via `loadamberprep` + `loadamberparams`, which is how we plug the chromophore parameters (published set or AM1-BCC) into the topology.
- AMBER's `cpptraj` provides the trajectory-analysis primitives (RMSD, RMSF, Rg, hbond, native contacts, SS retention, PDB extraction) needed for cross-candidate comparison.

> **Reminder:** Force fields are parameterized near room temperature. 352.15 K MD is a **comparative stress test**, not a measurement of T_m.

---

## 2. Required Input Files

### A. Structural inputs (one per system — 6 candidates + 1 reference WT)

| File | Source | Notes |
|---|---|---|
| `<candidate>_AF3.pdb` | AlphaFold 3 server | Single-chain GFP monomer with residues 65–67 still uncyclized |
| `<candidate>_clean.pdb` | `pdb4amber` on AF3 output | Clean atom names, no HETATM clutter, chromophore residues replaced by a single `CRO` residue (see §3) |
| `complex.pdb` (input filename used in `addH.in`) | derive from `<candidate>_clean.pdb` per system | The existing `addH.in` script reads a file literally named `complex.pdb`. Rename per system or symlink. **needs verification** that the chromophore residue in this PDB matches the library residue name used in step B. |

### B. Chromophore / nonstandard-residue parameter files

| File | Purpose | Status in this folder |
|---|---|---|
| `file.prepin` | Amber prep file for the cyclized GFP chromophore (`CRO`) | **Not present** in `Command/`. The existing `addH.in` and `addWater.in` reference `file.prepin` and `file.frcmod` as a generic ligand placeholder. For GFP, this must be the GFP chromophore prep from a **published Amber-compatible parameter set** (preferred), or generated with AM1-BCC charges via `antechamber -c bcc` (no QM engine required). **needs verification** — cite the parameter source and place alongside `complex.pdb`. |
| `file.frcmod` | Missing parameters for `CRO` (bonds/angles/torsions/vdW) | Same as above — must be the GFP chromophore frcmod, not a small-molecule frcmod. **needs verification**. |
| `atomic_ions.lib` | Ion library (referenced in `addALL.in`) | Standard AMBER file, supplied by AmberTools |
| `frcmod.ions234lm_1264_tip3p` | Ion parameters (referenced in `addALL.in`) | Standard AMBER file. **Suggested improvement:** for OPC water you may prefer the matching OPC ion set rather than tip3p ions; **needs verification** against your AmberTools version. |

### C. AMBER input files (already in `Command/`)

| File | Stage | Notes |
|---|---|---|
| `addH.in` | tleap, dry build | Loads `leaprc.protein.ff19SB`, `leaprc.gaff2`, ligand prep+frcmod, builds `complex_addH.{pdb,top,crd}`. **Contains disulfide bonds (`bond w.7.SG w.33.SG`, etc.) that do not apply to sfGFP — remove for GFP.** |
| `addWater.in` | tleap, solvation | Loads OPC water, adds 1 Na+ (literal count, not neutralization), solvates with `solvateoct OPCBOX 12.0`. **GFP adaptation: switch to `addIons w Na+ 0` + `addIons w Cl- 0` for charge neutralization, and add a physiological-ion step (see §4).** |
| `addALL.in` | alternative tleap script | More complete: uses `set default PBradii mbondi3`, `solvateoct OPCBOX 13 iso`, loads ion libs. Easier starting point for GFP than `addH.in`+`addWater.in`. **GFP adaptation: remove the three disulfide-bond lines (residues 24/35, 36/55, 46/89) — those do not exist in sfGFP.** |
| `min_H.in` | hydrogen-only minimization | Uses `ibelly=1` + an atom-by-atom `FIND` block over `RES 1 302`. **GFP: change `RES 1 302` → `RES 1 239`** (or whatever residue count your built system has after chromophore replacement — **needs verification**). |
| `minWAT.in` | minimize water/ions while restraining solute | `restraintmask` via belly `RES 1 302`. **GFP: change to `RES 1 <n_protein_residues>`** |
| `minALL.in` | full-system minimization (no restraints) | Already restraint-free; no residue-range change required. |
| `min1.i` | staged min, restrain non-water/ion heavy atoms | `restraintmask='!(:WAT,Na+,Cl-  | @H=)'` — OK for GFP. |
| `min2.i` | restrain backbone atoms (and ZN) | `restraintmask='!:WAT&@CA,C,O,N,ZN'` — **remove `ZN`** for GFP. |
| `min3.i` | restrain CA (and ZN) | `restraintmask='@CA,ZN'` — **remove `ZN`** for GFP. |
| `min4.i` | no restraints | OK for GFP as-is. |
| `heat.i` | NVT heat 10 K → 310 K over 200 ps | `restraintmask='!:WAT&@CA,N,C,O,ZN'` — **remove `ZN`**. For 79 °C work, additionally heat to 352.15 K (see §8). |
| `eq1.i` | NPT equilibration, 500 ps, restraint 5.0 kcal/mol/Å² | Same mask issue — **remove `ZN`**. |
| `eq2.i` | NPT equilibration, 500 ps, restraint 2.0 kcal/mol/Å² | Same mask issue — **remove `ZN`**. |
| `eq4.i` | 2 ns NPT equilibration, `ntr=0`, `cut=10.0`, **`barostat=1` (Berendsen)** | Present in folder but **not in the `min_heat.qsub` loop**. Optional extended equilibration step. **needs verification** whether you want to insert this between `eq2` and production. Note its `cut=10.0` differs from the rest of the workflow's `cut=12.0`. |
| `md.i` | 2 ns NPT production, no restraints, 310 K | Active production input. For 79 °C work see §8. |
| `hbond.in` | cpptraj hbond on `:1-302` | **GFP: change `:1-302` → `:1-239`** (or your protein range). |
| `mmpbsa.in` | MM-PBSA / MM-GBSA settings | Designed for receptor–ligand decomposition (`receptor_mask=":1-140"`, `ligand_mask=":141"`). **Does not directly apply to a single-chain GFP** unless you redefine a "ligand" (e.g. the chromophore residue) — discuss in §9. |

### D. Job/launch scripts (already in `Command/`)

| File | Stage | Notes |
|---|---|---|
| `minH.job` | runs `sander` for `min_H.in` | Uses serial `sander`. Suitable for the small belly minimization. |
| `minWAT.job` | runs `pmemd.MPI` (np=10) for `minWAT.in` | Adjust `-np` to your hardware. |
| `minALL.job` | runs `pmemd.MPI` (np=10) for `minALL.in` | Same. |
| `min_heat.qsub` | sequential `pmemd.cuda` loop over `min1 → min2 → min3 → min4 → heat → eq1 → eq2` | GPU required. Set `$PMEMDHOME` in your shell first. **Does not include `eq4.i`** — see notes above. |
| `runMD.qsub` | sequential `pmemd.cuda` loop, default `eq2 → md1 … md50` (≈100 ns); commented continuation blocks up to `md500` (≈1 µs) | GPU required. Loop variable `name` walks the input filenames. |
| `mmpbsa.job` | MMPBSA via `MMPBSA.py.MPI` | Tied to ligand–receptor masks; not directly used for GFP single-chain unless re-defined. |

### E. CPPTRAJ analysis (already in `Command/CPPTRAJ/`)

| File | Purpose | GFP-relevant notes |
|---|---|---|
| `cpptraj-analysis.sh` | master script — runs all 8 steps in order | Calls `hbond.in` which uses `:1-302`; update for GFP. |
| `sumtotal.in` | combine MD chunks, strip water/ions, write `sum_MD.nc` + `ref-complex.rst7` | Uses `complex_md{1..250}.nc`, `:1-141` for centering. **Adjust to your actual MD chunk count and to GFP residue range.** |
| `summmpbsa.in` | extract MMPBSA window `md{226..250}` | Adjust ranges to actual chunks; only relevant if doing MMPBSA. |
| `rmsd.in` | RMSD for complex/protein/ligand/backbone/binding-site | Uses `:1-302`, `:1-301`, `:302`. **For GFP single chain, drop "ligand" selection and replace with chromophore residue index.** |
| `rmsf.in` | RMSF per residue (backbone) + B-factor outputs | Residue-range-agnostic; works with `byres`. |
| `contact_5.in` | native contacts ligand `:302` ↔ protein `:1-301` | **For GFP single chain, redefine as chromophore (CRO residue index) ↔ rest of protein.** |
| `run-extract.in` | extract every 10th frame as PDB | Centers/aligns on `:1-140`. **Adjust to GFP protein range.** |
| `run-Last-PDB.in` | write last frame as PDB | Same centering note. |
| `sumMD-rms.in` | optional ASCII trajectory | Uses `:1-302`. Adjust. |

---

## 3. GFP Structure Preparation

### Input expected for AMBER

A clean, single-chain, hydrogen-friendly PDB with one chromophore residue (`CRO`) covering what was originally Thr65 / Tyr66 / Gly67. **AlphaFold 3 does NOT cyclize the chromophore** — it returns three standard residues at positions 65–67. You must manually build the cyclic chromophore before AMBER.

### Workflow

1. **Start from the AlphaFold 3 PDB** for the candidate (`<candidate>_AF3.pdb`).
2. **Inspect** in PyMOL or ChimeraX. Check:
   - β-barrel intact (11 strands + central helix carrying the chromophore).
   - No chain breaks or missing residues.
   - Per-residue pLDDT acceptable for the chromophore-microenvironment residues (Cα positions T65, Y66, G67, H148, T203, S205, R96, E222 using competition sfGFP-239 numbering).
3. **Build the cyclized chromophore.** Replace residues 65–66–67 with a single `CRO` residue. Two practical paths:
   - **(a) Use a published GFP chromophore PDB template** (e.g. PDB 2B3P or any sfGFP X-ray structure with `CRO`/`GYC`) and overlay residues 64–68 from the template onto your AlphaFold model. Save the merged model with consistent atom names matching the prep file you will use in `tleap`.
   - **(b) Generate parameters with AM1-BCC** via `antechamber -c bcc` on the extracted chromophore (no QM engine required), then `parmchk2 -s gaff2`. Heavier than (a) — only use if a matching published set is not available. **needs verification.**
4. **Decide protonation state of the phenolic OH on Tyr-derived ring** (A-form neutral vs B-form anionic). The B-form (deprotonated phenolate, total chromophore charge −1, multiplicity singlet) is the canonical bright state and is the recommended choice unless you have a specific reason to model the A-form. **Document the choice in writing** before running `tleap`. **needs verification** — the chromophore residue total charge in the prep file must match this choice.
5. **Run `pdb4amber`:**
   ```bash
   pdb4amber -i <candidate>_chromophorized.pdb -o <candidate>_clean.pdb
   ```
   This standardizes atom names, removes alt-locs, adds chain-end residues, flags any oddities in `pdb4amber.log`.
6. **Set HIS protonation** by inspection. AmberTools will default to HIE unless you rename (`HID`, `HIE`, `HIP`). H148 in particular is in the chromophore environment — choose deliberately. **needs verification** per candidate.
7. **No disulfide bonds for sfGFP.** sfGFP/avGFP have Cys residues but they do not form S–S bonds in the folded protein. **Remove** the `bond w.X.SG w.Y.SG` lines from `addH.in`/`addWater.in`/`addALL.in` for any GFP build.
8. **Save as `complex.pdb`** in the working directory (or change the `loadpdb` line in `addH.in` to your filename).

### Files to verify before running `tleap`

```bash
ls complex.pdb                # the GFP+CRO PDB
ls file.prepin                # chromophore prep (must be GFP CRO, not a generic ligand)
ls file.frcmod                # chromophore frcmod
ls addH.in addWater.in        # tleap scripts you intend to use
```

If `file.prepin` / `file.frcmod` are not the GFP chromophore parameters, **stop and obtain them first** — from a published Amber-compatible CRO parameter set (preferred) or via `antechamber -c bcc` + `parmchk2` (no QM needed). The names are placeholders carried over from the previous ligand-system workflow.

---

## 4. AMBER System Building (`tleap`)

### 4.1 Force-field choices fixed by the existing scripts

| Component | Choice | Source line |
|---|---|---|
| Protein force field | **ff19SB** | `source leaprc.protein.ff19SB` in `addH.in` / `addWater.in` / `addALL.in` |
| Chromophore force field for nonstandard residue | **GAFF2** | `source leaprc.gaff2` |
| Water model | **OPC** (4-point, recommended pairing with ff19SB) | `source leaprc.water.opc` in `addWater.in` / `addALL.in` |
| Ion library | `atomic_ions.lib` + `frcmod.ions234lm_1264_tip3p` | `addALL.in`. **Suggested improvement (needs verification):** switch to the OPC-matched ion frcmod if available in your AmberTools version. |
| PB radii (for downstream PBSA, if used) | `mbondi3` | `set default PBradii mbondi3` in `addWater.in` / `addALL.in` |

### 4.2 GFP-adapted `tleap` script (recommended starting point)

Take `addALL.in` as the cleaner template, then **strip ligand-system lines that don't apply to GFP**. A minimal GFP-ready script (do **not** overwrite the original — save as `addGFP.in`):

```bash
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.opc

# Chromophore parameters — GFP CRO prep+frcmod from a published set or antechamber AM1-BCC (no QM)
loadamberprep   ./file.prepin
loadamberparams ./file.frcmod

# Ion parameters
loadamberparams atomic_ions.lib
loadamberparams frcmod.ions234lm_1264_tip3p     # needs verification vs OPC-matched ions

set default PBradii mbondi3

w = loadpdb complex.pdb

# NO disulfide bonds in sfGFP — intentionally omitted

check w
charge w

# Neutralize charge (NOT a fixed-count addIons)
addIons w Na+ 0
addIons w Cl- 0

# Truncated octahedral OPC box, 12 Å buffer (matches the original workflow's cut=12.0)
solvateoct w OPCBOX 12.0

check w
charge w

savepdb w addWAT.pdb
saveamberparm w addWAT.top addWAT.crd

quit
```

### 4.3 Run and check

```bash
tleap -f addGFP.in 2>&1 | tee tleap.log
```

Then inspect `leap.log` / `tleap.log` and your saved files:

```bash
ls addWAT.pdb addWAT.top addWAT.crd
grep -Ei 'warning|error|fatal|missing|unknown' leap.log
```

Acceptance criteria:

- No `FATAL` / `Missing parameters` / `Unknown residue` lines.
- `check w` reports no impossibly close contacts or missing connectivity.
- `charge w` reports an integer total charge **before** the `addIons … 0` calls, and **0.0** after.
- `addWAT.pdb` opens cleanly in PyMOL/ChimeraX with a single chain, one CRO residue near the barrel center, and a truncated octahedron of waters.

> **GFP-specific note on the original `addWater.in`:** it uses `addIons w Na+ 1` — that adds exactly **one** sodium ion, regardless of net charge. For GFP at neutral pH the net protein charge is generally non-zero and a fixed `Na+ 1` will leave a non-neutral box. Use `addIons … 0` for neutralization (this is the recommended improvement).

---

## 5. Minimization Workflow

The existing workflow runs minimization in two scripted phases:

**Phase A — solvent + full-system minimization (CPU `pmemd.MPI` / `sander`)**
`addH → minH (belly H-only) → ambpdb → addWater → minWAT (solute frozen) → minALL (no restraints)`

**Phase B — staged restraint relaxation on GPU (`pmemd.cuda`), inside `min_heat.qsub`**
`min1 → min2 → min3 → min4`

### 5.1 Phase A — Solvent / full-system minimization

#### Hydrogen-only minimization

| Item | Value |
|---|---|
| Input | `min_H.in` |
| Job | `minH.job` (calls `sander`) |
| Command | `bash minH.job` |
| Reads | `complex_addH.top`, `complex_addH.crd` |
| Writes | `complex_minH.restrt`, `min_H.out` |
| Bridge | `ambpdb -p complex_addH.top -c complex_minH.restrt > minH.pdb` |

**GFP adaptation:** change `RES 1 302` → `RES 1 <n_residues_in_complex_addH.pdb>` in `min_H.in`. **needs verification** that all `FIND ... SEARCH` H-name groups in `min_H.in` cover the GFP atom-name set (in particular `HA`, `HB`, `HG`, `H1`–`H5`). The block in the existing file is exhaustive and should be fine.

#### Solvent (water + ions) minimization, solute restrained

| Item | Value |
|---|---|
| Input | `minWAT.in` |
| Job | `minWAT.job` (`mpirun -np 10 pmemd.MPI`) |
| Command | `bash minWAT.job` |
| Reads | `addWAT.top`, `addWAT.crd` |
| Writes | `minWAT.restrt`, `minWAT.out` |

**GFP adaptation:** in `minWAT.in`, change `RES 1 302` → `RES 1 <n_protein_residues>`.

#### Full-system minimization (no restraints)

| Item | Value |
|---|---|
| Input | `minALL.in` |
| Job | `minALL.job` |
| Command | `bash minALL.job` |
| Reads | `addWAT.top`, `minWAT.restrt` |
| Writes | `min_all.restrt`, `minALL.out` |

No GFP-specific edits needed (already `ntr=0`).

#### Quality checks after Phase A

```bash
tail -n 40 min_H.out minWAT.out minALL.out
ls minH.pdb complex_minH.restrt minWAT.restrt min_all.restrt
```

Look for: monotone energy decrease, no `NaN`, no SHAKE warnings (none expected for `imin=1`), reasonable final RMS gradient.

### 5.2 Phase B — Staged restraint relaxation (GPU)

Driven by `min_heat.qsub`, which copies `min_all.restrt → complex_ini.rst7` then walks the loop `for name in min1 min2 min3 min4 heat eq1 eq2`.

| Stage | Input file | Restraints | Purpose |
|---|---|---|---|
| `min1` | `min1.i` | `!(:WAT,Na+,Cl- | @H=)`, 25 kcal/mol/Å² | Free water, ions, and H; freeze the rest |
| `min2` | `min2.i` | `!:WAT&@CA,C,O,N,ZN`, 25 kcal/mol/Å² | Relax side chains; keep backbone (and ZN — **remove for GFP**) restrained |
| `min3` | `min3.i` | `@CA,ZN`, 25 kcal/mol/Å² | Restrain only Cα (remove ZN) |
| `min4` | `min4.i` | none | Free full-system minimization |

**Required GFP edits to `Command/*.i`:**

- `min2.i`, `min3.i`, `heat.i`, `eq1.i`, `eq2.i`: **remove `,ZN`** from every `restraintmask`. GFP has no zinc.

After Phase B you should have:

```
complex_min1.rst7 complex_min2.rst7 complex_min3.rst7 complex_min4.rst7
complex_min1.mdout ... complex_min4.mdout
```

Spot-check:

```bash
for n in min1 min2 min3 min4; do
  echo "=== $n ==="; tail -n 30 complex_${n}.mdout
done
```

Look for: clean termination, monotone energy decrease, no `NaN`, no `vlimit exceeded`.

---

## 6. Heating and Equilibration

These stages also run inside `min_heat.qsub`, immediately after `min4`.

### 6.1 Heating

| Item | Value |
|---|---|
| Input | `heat.i` |
| Ensemble | **NVT** (`ntb=1, ntp=0`) |
| Thermostat | Langevin, `ntt=3`, `gamma_ln=5.0` |
| Time | 100 000 × 2 fs = **200 ps** |
| Temperature ramp | 10 K → **310 K** (via `nmropt=1` + `&wt TEMP0`) |
| Restraints | `restraintmask='!:WAT&@CA,N,C,O,ZN'`, 5.0 kcal/mol/Å² |
| SHAKE | on (`ntc=2, ntf=2`) |
| Cutoff | 12.0 Å (consistent with OPC) |

**GFP edits:** remove `,ZN`. If you want the production at 352.15 K (see §8), keep this 300-K heat as-is and add a second heating block 310 K → 352.15 K before the high-T production run.

### 6.2 Equilibration 1 — NPT, strong backbone restraint

| Item | Value |
|---|---|
| Input | `eq1.i` |
| Ensemble | **NPT** (`ntb=2, ntp=1`) |
| Barostat | Monte Carlo, `barostat=2`, `taup=2.0`, `pres0=1.01325` |
| Thermostat | Langevin, `gamma_ln=2.0`, `temp0=310.0` |
| Time | 250 000 × 2 fs = **500 ps** |
| Restraints | `!:WAT&@CA,N,C,O,ZN`, **5.0** kcal/mol/Å² |

**GFP edits:** remove `,ZN`.

### 6.3 Equilibration 2 — NPT, reduced restraint

| Item | Value |
|---|---|
| Input | `eq2.i` |
| Ensemble | NPT |
| Time | 500 ps |
| Restraints | same mask, **2.0** kcal/mol/Å² |

**GFP edits:** remove `,ZN`.

### 6.4 Optional `eq4.i`

`eq4.i` is **present in the folder but not used by `min_heat.qsub`**. It is a 2 ns NPT MD with `ntr=0`, `cut=10.0`, `barostat=1` (Berendsen). If you want an additional unrestrained equilibration before production, insert `eq4` into the loop:

```bash
# Suggested improvement — needs verification before use
for name in min1 min2 min3 min4 heat eq1 eq2 eq4 ; do
```

**needs verification:** `eq4.i` uses `cut=10.0`, which differs from the workflow's `cut=12.0`. Either unify to 12.0 or document the change.

### 6.5 Stability checks after heating + equilibration

```bash
for n in heat eq1 eq2; do
  echo "=== $n ==="
  tail -n 50 complex_${n}.mdout
done
```

Look for:

- Temperature reaches and holds 310 K (heat: rises smoothly; eq1/eq2: fluctuates around 310 ± a few K).
- Pressure: fluctuates in NPT but does not diverge; mean ≈ 1 bar.
- Density: settles near 1.00 g/cm³ (slightly lower for OPC at 310 K).
- Total energy: no sudden non-physical jumps.
- No `SHAKE`, `vlimit`, or `NaN` messages.

If any check fails, **do not** proceed to production. Re-minimize / re-equilibrate with stronger initial restraints or a shorter time-step.

---

## 7. Production MD Simulation

### 7.1 Existing production setup

| Item | Value |
|---|---|
| Input | `md.i` |
| Launcher | `runMD.qsub` |
| Engine | `pmemd.cuda` |
| Ensemble | NPT |
| Temperature | `temp0=310.0` |
| Time per chunk | 1 000 000 × 2 fs = **2 ns** |
| Active default loop | `eq2 → md1 → md2 → … → md50` (≈100 ns total) |
| Commented continuation blocks | `md51 → md250` (extra 400 ns) and `md251 → md500` (extra 500 ns), for up to 1 µs |
| Restart chain | `old=$name` after each chunk — each chunk reads the previous `.rst7` |

### 7.2 How to run

```bash
# Pre-flight
ls addWAT.top complex_eq2.rst7 min_all.restrt md.i runMD.qsub

# Local
bash runMD.qsub

# Or on a scheduler (PBS / Torque)
qsub runMD.qsub
```

### 7.3 How restart chaining works (re-stated for clarity)

- First iteration: `old=eq2`, `name=md1` → reads `complex_eq2.rst7`, writes `complex_md1.{nc,rst7,mdout,mdinfo}`.
- Loop sets `old=md1` and proceeds to `md2`, and so on.
- If any chunk fails, its `.rst7` is not produced, and subsequent chunks cannot start.

### 7.4 How to continue past `md50`

Uncomment one of the continuation blocks in `runMD.qsub` and update `old`. Cleaner alternative using `seq` (already suggested in `AMBER-Guide.md`):

```bash
start=51 ; end=250 ; old=md50
for i in $(seq $start $end); do
  name=md${i}
  $PMEMDHOME/bin/pmemd.cuda -O -i md.i \
    -o complex_${name}.mdout -p addWAT.top -c complex_${old}.rst7 \
    -ref min_all.restrt -x complex_${name}.nc -r complex_${name}.rst7 \
    -inf complex_${name}.mdinfo
  old=${name}
done
```

### 7.5 Per-chunk output

| File | Meaning |
|---|---|
| `complex_mdX.nc` | NetCDF trajectory chunk (2 ns at `ntwx=5000`) |
| `complex_mdX.rst7` | Restart for chunk `X+1` |
| `complex_mdX.mdout` | Energies / temperature / pressure log |
| `complex_mdX.mdinfo` | Live runtime info |

### 7.6 Monitoring

```bash
tail -f complex_md1.mdinfo            # live status
tail -n 50 complex_md1.mdout          # post-chunk summary
ls complex_md*.nc | wc -l             # chunks completed
```

---

## 8. Thermostability Testing (79 °C / 352.15 K)

### 8.1 Temperature conversion

```
T(K) = T(°C) + 273.15
79 °C  → 352.15 K
72 °C  → 345.15 K   (the SynBio 2026 assay temperature per ml_part/README.md)
```

Pick one as primary, document it in the methods section, and apply it to **all** systems identically.

### 8.2 Two practical options (do not invent new files; copy and edit existing ones)

#### Option A — Sequential heat-up then high-T production *(recommended)*

1. Run the existing protocol up to `complex_eq2.rst7` at 310 K (Phase A → Phase B → heat → eq1 → eq2).
2. **Copy `heat.i` to `heat_hi.i`** and change the `&wt` ramp to 310 → 352.15:
   ```
   &wt
     TYPE='TEMP0', ISTEP1=0, ISTEP2=100000,
     VALUE1=310.0, VALUE2=352.15,
   /
   &wt TYPE='END' /
   ```
   Also set `tempi=310.0` and keep `nmropt=1`.
3. **Copy `md.i` to `md_352.i`** and change:
   ```
   tempi=352.15,
   temp0=352.15,
   ```
4. Modify `runMD.qsub` (or make `runMD_352.qsub`) to run `heat_hi` once, then loop `md1_352 → md2_352 → …`.
5. Run at least **2 independent replicates per system** with different `ig` seeds (`ig=-1` already randomizes, but record the chosen value from `mdout` for traceability).

#### Option B — Direct high-T production after extended 310 K equilibration

Use the existing `md.i` but flip both `tempi` and `temp0` to 352.15. **Riskier** — the thermostat will have to absorb a 42 K jump at chunk start. Not recommended for production-quality data.

### 8.3 Reference / control system

Run **WT sfGFP-239** through the identical pipeline (same chromophore parameters, same buffer thickness, same ion count, same heat-up, same chunk count, same replicates, same seeds). Without this control you cannot tell whether observed instability comes from the candidate's design or from your simulation setup.

### 8.4 Fair-comparison rules

- Identical force field, water model, ion model.
- Identical box buffer (e.g. `solvateoct OPCBOX 12.0`).
- Identical minimization → heat → equilibration → high-T heat protocol.
- Identical production length and replicate count.
- Identical CPPTRAJ analysis settings and atom selections (see §9).
- WT control included in every plot and table.

### 8.5 Wording in results

Use cautious language. Examples:

- ✅ *"Candidate 3 shows lower mean Cα-RMSD than WT sfGFP under the 50 ns 352.15 K NPT protocol with ff19SB / OPC."*
- ✅ *"This is consistent with relatively greater structural stability in silico under these conditions."*
- ✅ *"Experimental measurement (e.g. F_final/F_initial after the assay heat challenge) is required to confirm this prediction."*
- ❌ *"Candidate 3 is thermostable at 79 °C."*
- ❌ *"MD proves the design will work."*

---

## 9. CPPTRAJ Analysis

The master script `Command/CPPTRAJ/cpptraj-analysis.sh` runs 8 steps in order. **GFP adaptation: residue ranges and the "ligand" definition need to be edited** because the existing files were written for a 302-residue protein–ligand complex.

### 9.1 Step-by-step pipeline (with GFP edits flagged)

#### Step 1 — Combine production chunks and strip solvent (`sumtotal.in`)

- Reads `addWAT.top` and `complex_md{1..250}.nc`.
- `strip :WAT,Na+,Cl-`, `autoimage`, writes `sum_MD.nc` and `stripped.addWAT.top`.
- Writes `ref-complex.rst7` from `min_all.restrt`.

**GFP edits:**
- Update the trajectory range to your actual chunk count: `trajin complex_md{1..50}.nc` (for 100 ns) or `{1..N}` for whatever you ran.
- Change centering selection `:1-141` → **`:1-<n_protein_residues>`** for GFP (e.g. `:1-238` for sfGFP-239 where CRO is residue 239, or `:1-237` if CRO replaces a triad in the middle — **needs verification** based on your final residue numbering after chromophore replacement).
- The `surf` selection of "binding-site" residues is **not applicable** to GFP — remove it or replace with a GFP-relevant selection (e.g. chromophore environment H148/T203/S205/R96/E222 in your built residue numbering).

#### Step 2 — Prepare MMPBSA window (`summmpbsa.in`)

Only relevant for **interaction-energy** analysis (e.g. chromophore ↔ rest-of-protein, if you want to decompose). For pure thermostability assessment, **skip this step**.

#### Step 3 — RMSD (`rmsd.in`)

The existing file calculates 5 RMSDs (complex / protein / ligand / backbone / binding-site). For GFP, replace with:

```
parm stripped.addWAT.top
trajin sum_MD.nc

center :1-<NPROT> mass origin
image origin center
autoimage familiar

reference ref-complex.rst7

# Backbone RMSD over the whole monomer
rms reference :1-<NTOTAL>@CA,C,O,N out rmsd_backbone.dat
# β-barrel Cα only (provide your barrel residue list — needs verification per AF model)
rms reference :<BARREL_LIST>@CA   out rmsd_barrel.dat
# Chromophore-region heavy atoms (residue index of CRO in your topology — needs verification)
rms reference :<CRO_IDX>&!@H=     out rmsd_cro.dat
```

Where `<NPROT>` is the protein residue count, `<NTOTAL>` includes CRO, `<BARREL_LIST>` is the explicit Cα selection of the 11 β-strands (define once, reuse across all systems), and `<CRO_IDX>` is the residue index of your `CRO` after build (**needs verification per system — print with `resinfo` in cpptraj after parm load**).

#### Step 4 — RMSF (`rmsf.in`)

Per-residue backbone RMSF + B-factor outputs. Works for GFP without selection changes — just confirm the system residue count matches and use the same reference atoms across all systems for fairness.

#### Step 5 — Hydrogen bonds (`hbond.in`)

Existing file uses `:1-302`. **GFP edit:**

```
parm stripped.addWAT.top
trajin sum_MD.nc
autoimage
reference ref-complex.rst7

hbond :1-<NTOTAL> nointramol dist 3.5 angle 120.0 \
  out nhb_all.dat avgout avghb_all.dat
```

Adjust `<NTOTAL>` and the `nointramol` flag depending on whether you want only intra-protein or include the chromophore-to-protein H-bonds (chromophore H-bond network is informative for GFP stability).

#### Step 6 — Native contacts (`contact_5.in`)

Existing file does ligand `:302` ↔ protein `:1-301`. For GFP, useful re-definitions:

- **(a) Chromophore ↔ rest of protein:** `nativecontacts :<CRO_IDX>&!@H= :1-<NPROT>&!@H= … distance 5.0`
- **(b) β-barrel internal native contacts:** `nativecontacts :<BARREL_LIST>&@CA :<BARREL_LIST>&@CA … distance 8.0` to track barrel integrity over the trajectory.

#### Step 7 — Radius of gyration *(recommended addition, not in folder)*

Suggested improvement — create `radgyr.in`:

```
parm stripped.addWAT.top
trajin sum_MD.nc
autoimage
radgyr :1-<NPROT>@CA out rg.dat
```

Increasing Rg is a classical unfolding indicator.

#### Step 8 — Secondary structure *(recommended addition)*

Suggested improvement — create `secstruct.in`:

```
parm stripped.addWAT.top
trajin sum_MD.nc
autoimage
secstruct :1-<NPROT> out ss.dat sumout ss_sum.dat
```

Loss of β-strand content at 352.15 K is a strong "barrel opening" indicator.

#### Step 9 — Extract representative PDBs

Use the existing `run-extract.in` (every 10th frame) and `run-Last-PDB.in` (last frame). **GFP edit:** change `:1-140` → `:1-<NPROT>` for centering/alignment.

#### Step 10 — Optional MM/PBSA

`mmpbsa.in` + `mmpbsa.job` are written for a receptor (`:1-140`) + ligand (`:141`) decomposition. For a GFP single chain, MMPBSA does not have an obvious interpretation **unless** you redefine the "ligand" as the chromophore residue and the "receptor" as the rest of the protein, to estimate a chromophore–protein interaction energy. **Treat as exploratory; needs verification.**

### 9.2 Running the analysis

```bash
cd Command/CPPTRAJ
bash cpptraj-analysis.sh
# inspect each .log for errors
for f in sumtotal summmpbsa rmsd rmsf hbond contact_5 run-extract run-Last-PDB; do
  echo "=== $f ==="; tail -n 30 ${f}.log
done
```

### 9.3 Cross-candidate output table

Build a single combined CSV per stage (300 K, 352.15 K) including the WT row:

```
candidate,replicate,mean_BB_RMSD_A,max_BB_RMSD_A,mean_Rg_A,frac_HB_retained,frac_strand_retained,mean_CRO_env_dist_A
ref_sfGFP,1,...
cand_01,1,...
cand_01,2,...
...
```

Only enter values you actually computed — never fabricate numbers.

---

## 10. Quality Control and Troubleshooting

### 10.1 `tleap` errors

| Symptom | Likely cause | Fix |
|---|---|---|
| `Unknown residue: CRO` (or similar chromophore name) | `loadamberprep file.prepin` not run, or prep file uses a different residue name than the PDB | Confirm residue name matches between `complex.pdb` and the prep file; reload prep before `loadpdb` |
| `Missing parameters` / `*_*` torsion not found | `file.frcmod` incomplete | Re-run `parmchk2 -s gaff2 -i chromophore.mol2 -o file.frcmod`; supplement missing terms from a published GFP parameter set and cite it |
| `FATAL: Could not find atom XXXX` | Atom-name mismatch between PDB and prep | Standardize names; re-run `pdb4amber`; manually rename if needed |
| Final `charge w` not integer | Wrong protonation state on a residue or wrong total chromophore charge | Inspect HID/HIE/HIP assignments; re-check chromophore A/B-form choice |
| Non-zero charge after `addIons … 0` calls | Library/prep charges inconsistent | Recheck chromophore charges (published set / re-run `antechamber -c bcc`); rebuild |

### 10.2 Minimization / heating errors

| Symptom | Likely cause | Fix |
|---|---|---|
| `SHAKE failure` early in `heat.i` | Bad starting geometry, clashes, or too-fast heating | Run additional minimization; temporarily reduce `dt` to 0.001; lengthen the temperature ramp |
| `vlimit exceeded` | Numerical instability — usually clash | Inspect last frame visually; rebuild |
| Energy `NaN` in early MD | Almost always a build issue (missing parameters, mis-named atom, bad chromophore connectivity) | Go back to `tleap` and `pdb4amber`; do **not** "fix" by running more MD |
| Density not converging | Insufficient equilibration / box too small | Extend `eq1`/`eq2`; consider a larger buffer (try `OPCBOX 13`) |
| Temperature overshoot during heat | Restraint too tight or barostat conflicting with NVT | Verify `ntb=1, ntp=0` in heat.i; check `gamma_ln` |

### 10.3 Production / restart-chain errors

| Symptom | Likely cause | Fix |
|---|---|---|
| `complex_eq2.rst7` not found | Equilibration didn't finish | Inspect `complex_eq2.mdout`; rerun `min_heat.qsub` from the last good restart |
| `complex_mdN.rst7` missing mid-chain | Chunk N crashed | Check `complex_mdN.mdout`; restart that chunk from `complex_md(N-1).rst7` |
| Trajectory chunks have wrong number of atoms in `cpptraj` | Topology and trajectory mismatch (mixed `addWAT.top` vs `stripped.addWAT.top`) | Use `parm` and `trajin` from the **same** build; never reload a different topology mid-analysis |

### 10.4 GFP-specific gotchas

| Issue | Why it happens | What to do |
|---|---|---|
| `bond w.X.SG w.Y.SG` failure in tleap | The existing scripts have placeholder disulfide bonds for the original ligand system | Remove those lines for any sfGFP build |
| `ZN` in restraint masks | Carried over from the original metal-binding system | Delete `,ZN` from `min2.i`, `min3.i`, `heat.i`, `eq1.i`, `eq2.i` |
| Chromophore "explodes" in heat.i | The chromophore prep/frcmod is wrong (often: total charge mismatch with the protonation state, or missing imidazolinone-ring torsions) | Stop and rebuild chromophore parameters from the published set / `antechamber`+`parmchk2`; do not bandage |
| AlphaFold model has 3 residues at positions 65–67 instead of `CRO` | AF3 does not cyclize the chromophore | Insert the cyclic chromophore from a published template before tleap (see §3) |
| Net charge non-zero after `addIons w Na+ 1` | The original script adds exactly 1 Na+ rather than neutralizing | Switch to `addIons w Na+ 0` + `addIons w Cl- 0` |

### 10.5 CPPTRAJ mask errors

| Symptom | Likely cause | Fix |
|---|---|---|
| `Atom selection produced 0 atoms` | Residue range or atom name wrong for your topology | `parm stripped.addWAT.top; resinfo; atominfo :<n>` inside cpptraj to inspect |
| RMSD trace is constant zero | Reference and trajectory loaded with different parms or different stripping | Use the `ref-complex.rst7` produced by the matching `sumtotal.in` run |
| `surf` failing on missing residues | Selection inherited from ligand-system version (e.g. residue 277, 282) which doesn't exist in GFP | Remove or replace the residue list |

### 10.6 Overinterpretation guards

- Never present an MD-derived ranking as definitive thermostability.
- Always show simulation length, replicate count, force-field versions, software versions.
- Always include the WT reference in every comparison plot / table.
- Use multi-metric rankings (RMSD + Rg + Hbond + SS retention + chromophore environment), not RMSD alone.

---

## 11. Final Workflow Summary

### 11.1 Complete command order (per candidate)

```bash
# --- Stage 1: dry build ---
tleap -f addGFP.in            # GFP-adapted version of addH.in / addALL.in

# --- Stage 2: hydrogen minimization (only if you split build from solvation) ---
bash minH.job
ambpdb -p complex_addH.top -c complex_minH.restrt > minH.pdb

# --- Stage 3: solvate, then solvent + full minimization ---
tleap -f addWater.in          # GFP-adapted (no disulfides; addIons … 0)
bash minWAT.job
bash minALL.job

# --- Stage 4: staged GPU minimization + heating + equilibration ---
bash min_heat.qsub            # min1 → min2 → min3 → min4 → heat → eq1 → eq2
                              # (optional: append eq4 if desired)

# --- Stage 5a: 300 K production (baseline) ---
bash runMD.qsub               # eq2 → md1 … md50  (≈100 ns; longer with continuation blocks)

# --- Stage 5b: high-T production at 352.15 K ---
bash heat_hi.qsub             # heat 310 → 352.15 K  (suggested improvement; needs verification)
bash runMD_352.qsub           # produce 50 ns–100 ns at 352.15 K, ≥2 replicates

# --- Stage 6: analysis ---
cd Command/CPPTRAJ
bash cpptraj-analysis.sh      # GFP-adapted masks
```

### 11.2 File-dependency map

```
addGFP.in / addH.in / addWater.in
  │
  ▼  tleap
addWAT.top, addWAT.crd, addWAT.pdb
  │
  ▼  minWAT.job → minALL.job
min_all.restrt
  │
  ▼  min_heat.qsub
complex_min1.rst7 → complex_min2.rst7 → complex_min3.rst7 → complex_min4.rst7
  │
  ▼
complex_heat.rst7 → complex_eq1.rst7 → complex_eq2.rst7
  │
  ▼  runMD.qsub
complex_md1.{nc,rst7,mdout,mdinfo} → ... → complex_mdN.{nc,rst7,mdout,mdinfo}
  │
  ▼  cpptraj-analysis.sh
sum_MD.nc, ref-complex.rst7, stripped.addWAT.top
  │
  ▼  rmsd.in, rmsf.in, hbond.in, contact_5.in, run-extract.in, run-Last-PDB.in
*.dat metrics + All-PDB/PDB.pdb.* + Last-NS.pdb
```

### 11.3 Pre-MD checklist

- [ ] AlphaFold model passed quality review (β-barrel intact, chromophore region acceptable).
- [ ] Chromophore residue built (cyclized; A/B-form chosen and documented).
- [ ] `file.prepin` and `file.frcmod` are the GFP chromophore parameters, not generic ligand placeholders.
- [ ] `complex.pdb` opens cleanly; `pdb4amber.log` clean.
- [ ] HIS residues set to HID/HIE/HIP deliberately.
- [ ] All disulfide-bond lines removed from `addH.in` / `addWater.in` / `addGFP.in` for sfGFP.
- [ ] All `ZN` strings removed from `min2.i`, `min3.i`, `heat.i`, `eq1.i`, `eq2.i`.
- [ ] `RES 1 302` in `min_H.in` and `minWAT.in` changed to GFP residue count.
- [ ] `addIons … 0` neutralization used (not fixed `Na+ 1`).
- [ ] Force fields/water/ions identical across all 7 systems (6 candidates + WT reference).
- [ ] `tleap.log` clean, total charge zero after ions.

### 11.4 Per-stage post-run checklist

**After minimization**
- [ ] Energies decrease, no `NaN`, no SHAKE messages.
- [ ] `min_all.restrt` and `complex_min{1..4}.rst7` exist.

**After heating**
- [ ] Temperature reaches 310 K smoothly.
- [ ] No `SHAKE` / `vlimit` / `NaN`.

**After equilibration**
- [ ] Temperature, pressure, density stable in eq2.
- [ ] `complex_eq2.rst7` exists.
- [ ] (If using eq4) extended NPT also stable.

**After production**
- [ ] All expected `complex_mdX.{nc,rst7,mdout}` exist.
- [ ] No crashes mid-chain; replicate trajectories distinguishable by file name.
- [ ] Same chunk count across all systems (WT included).

**After high-T production (if running 352.15 K)**
- [ ] Heat-up to 352.15 K gradual (not a 42 K jump).
- [ ] Same chunk count across all systems.
- [ ] At least 2 independent replicates per system, with seeds recorded.

### 11.5 Pre-analysis checklist

- [ ] Same atom-selection masks applied across all systems.
- [ ] Reference WT included in every comparison.
- [ ] Trajectory ranges in `sumtotal.in` / `summmpbsa.in` match actual chunk counts.
- [ ] Stripped topology and stripped trajectory came from the **same** build.
- [ ] Ranking uses multiple metrics, not RMSD alone.
- [ ] Limitations and "needs experimental validation" statement included in every output table.

---

## Distinguishing "existing" from "suggested improvement"

| Item | Existing in `Command/` | GFP-specific edit needed | Suggested improvement (not in folder) |
|---|---|---|---|
| `addH.in` build | yes | remove disulfide bonds; rename ligand prep to GFP CRO prep | use `addALL.in` instead — more complete |
| `addWater.in` | yes | remove disulfide bonds; switch `Na+ 1` → `Na+ 0` + `Cl- 0` | add physiological 0.15 M ions |
| Disulfide bonds in tleap | yes (residues 7/33, 19/130, 77/87 in `addH.in`; 24/35, 36/55, 46/89 in `addALL.in`) | **remove all** for sfGFP | — |
| `ZN` in restraint masks | yes (`min2.i`, `min3.i`, `heat.i`, `eq1.i`, `eq2.i`) | **remove all** | — |
| `RES 1 302` belly groups | yes (`min_H.in`, `minWAT.in`) | change to GFP residue count | — |
| `min1.i`–`min4.i` staged minimization | yes | none beyond ZN | — |
| `heat.i` 10 K → 310 K | yes | remove ZN | add `heat_hi.i` for 310 → 352.15 K |
| `eq1.i`, `eq2.i` | yes | remove ZN | optionally add `eq4.i` to the loop in `min_heat.qsub` |
| `md.i` production at 310 K | yes | none | add `md_352.i` for 352.15 K |
| `runMD.qsub` | yes | replace with `seq`-based loop for readability (optional) | add `runMD_352.qsub` for high-T run |
| `cpptraj-analysis.sh` 8 steps | yes | update masks (`:1-302` → GFP range), drop ligand-specific selections | add `radgyr.in` and `secstruct.in` |
| MMPBSA scripts | yes | inapplicable as-is to single-chain GFP | optionally re-cast chromophore↔protein decomposition |

---

## Next Actions

1. Obtain the GFP chromophore (`CRO`) prep + frcmod files from a **published Amber-compatible parameter set** (preferred) or via `antechamber -c bcc` + `parmchk2` (no QM engine). Place as `file.prepin` and `file.frcmod` in the working directory, and cite the source. **needs verification.**
2. Create the GFP-adapted versions of the scripts: `addGFP.in`, `min_H.in`, `minWAT.in`, `min2.i`, `min3.i`, `heat.i`, `eq1.i`, `eq2.i`, and (for thermostability) `heat_hi.i`, `md_352.i`, `runMD_352.qsub`. **Do not overwrite the originals** — save with new names.
3. Test the full pipeline on the **WT sfGFP-239** reference first. If WT builds and runs cleanly, the 6 candidates will follow the same pipeline with high confidence.
4. Decide and document the **primary thermostability temperature** (79 °C / 352.15 K vs 72 °C / 345.15 K) in the project methods.
5. After WT smoke test, run all 7 systems with ≥ 2 replicates per system, in two temperature blocks (300 K and high-T), keeping every protocol setting identical across systems.
6. Apply the CPPTRAJ pipeline with GFP-adapted masks; build the cross-candidate stability table; rank using multiple metrics, with WT in every comparison; report with cautious language and an explicit experimental-validation caveat.
