# AMBER Guide Following My Running Command Workflow

**Workflow focus:**  
This guide follows my AMBER preparation, minimization, heating, equilibration, production, and analysis workflow.

The preparation order is:

```text
tleap-H → minH → ambpdb → tleap Water → minWAT → minALL → min1 → min2 → min3 → min4 → heat → eq1 → eq2 → production MD → cpptraj analysis
```

This guide is organized around my actual file names and command style:

```text
addH.in
min_H.in
minH.job
addWater.in
minWAT.in
minWAT.job
minALL.in
minALL.job
min1.i
min2.i
min3.i
min4.i
heat.i
eq1.i
eq2.i
min_heat.qsub
runMD.qsub
cpptraj-analysis.sh
```

> **Scientific caution:** AMBER MD simulations can support interpretation of protein-ligand stability, conformational dynamics, and interaction persistence. However, MD simulation alone does not prove biological activity, inhibition, or therapeutic efficacy. Experimental validation is required for biological claims.

---

## 1. Complete Workflow Overview

The full workflow is divided into four main stages:

```text
Stage 1: Dry complex preparation
complex.pdb + file.prepin + file.frcmod
        ↓ addH.in / tleap
complex_addH.top + complex_addH.crd + complex_addH.pdb

Stage 2: Hydrogen minimization
complex_addH.top + complex_addH.crd
        ↓ min_H.in / minH.job
complex_minH.restrt
        ↓ ambpdb
minH.pdb

Stage 3: Solvation and early minimization
minH.pdb
        ↓ addWater.in / tleap
addWAT.top + addWAT.crd + addWAT.pdb
        ↓ minWAT.in / minWAT.job
minWAT.restrt
        ↓ minALL.in / minALL.job
min_all.restrt

Stage 4: Staged minimization, heating, and equilibration
min_all.restrt
        ↓ min_heat.qsub
complex_min1.rst7
complex_min2.rst7
complex_min3.rst7
complex_min4.rst7
complex_heat.rst7
complex_eq1.rst7
complex_eq2.rst7

Stage 5: Production MD
complex_eq2.rst7
        ↓ runMD.qsub
complex_md1.nc, complex_md1.rst7
complex_md2.nc, complex_md2.rst7
...
```

---

## 2. Required Input Files

Before running the workflow, check that these files exist:

```bash
ls complex.pdb
ls file.prepin
ls file.frcmod

ls addH.in
ls min_H.in
ls minH.job

ls addWater.in
ls minWAT.in
ls minWAT.job
ls minALL.in
ls minALL.job

ls min1.i
ls min2.i
ls min3.i
ls min4.i
ls heat.i
ls eq1.i
ls eq2.i
ls min_heat.qsub
```

---

## 3. Stage 1 — `tleap-H`: Add Hydrogens and Build Dry Complex

### Purpose

This stage creates a dry protein-ligand complex topology and coordinate file before water is added.

It uses:

- `leaprc.protein.ff19SB`
- `leaprc.gaff2`
- `file.prepin`
- `file.frcmod`
- `complex.pdb`
- manually defined disulfide bonds

### Input: `addH.in`

```bash
#!/bin/bash

source leaprc.protein.ff19SB
source leaprc.gaff2

loadamberprep   file.prepin
loadamberparams file.frcmod

w = loadpdb complex.pdb

# --- Disulfide bonds ---
bond w.7.SG w.33.SG
bond w.19.SG w.130.SG
bond w.77.SG w.87.SG

check w
charge w

savepdb w complex_addH.pdb
saveamberparm w complex_addH.top complex_addH.crd

quit
```

### Run

```bash
tleap -f addH.in
```

### Outputs

```text
complex_addH.pdb
complex_addH.top
complex_addH.crd
```

### Quality Check

```bash
ls complex_addH.pdb complex_addH.top complex_addH.crd
less leap.log
```

Check for:

- unknown residues
- missing ligand parameters
- incorrect ligand residue name
- failed disulfide bonds
- unusual charge
- missing atoms

---

## 4. Stage 2 — `minH`: Hydrogen Minimization

### Purpose

This step minimizes hydrogen atoms before solvation.

This is useful because hydrogens are often newly added by `tleap` or structure-preparation tools and may have strained positions.

### Input: `min_H.in`

```text
steps of minimization only H atoms
 &cntrl
 imin=1,
 maxcyc=3000,
 ncyc=1000,
 cut=12.0,
 ntpr=5,
 ntb=0,
 ibelly=1,
 &end
group input:select H
...
END
END
```

### Important Syntax Note

Your original `min_H.in` should include a comma after `ncyc=1000` if it is written on the same namelist block. A safer corrected namelist is:

```text
 &cntrl
 imin=1,
 maxcyc=3000,
 ncyc=1000,
 cut=12.0,
 ntpr=5,
 ntb=0,
 ibelly=1,
 &end
```

### Run: `minH.job`

```bash
#!/bin/bash

sander -O \
  -i min_H.in \
  -o min_H.out \
  -p complex_addH.top \
  -c complex_addH.crd \
  -r complex_minH.restrt \
  -ref complex_addH.crd
```

Run:

```bash
bash minH.job
```

### Output

```text
complex_minH.restrt
min_H.out
```

### Convert Restart to PDB

This step is essential because `addWater.in` reads `minH.pdb`.

```bash
ambpdb -p complex_addH.top -c complex_minH.restrt > minH.pdb
```

Check:

```bash
ls minH.pdb
```

---

## 5. Stage 3 — `tleap Water`: Add Water and Ions

### Purpose

This stage reloads the hydrogen-minimized dry complex and creates the solvated system.

It uses:

- `ff19SB`
- `GAFF2`
- `OPC` water model
- `mbondi3` radii for later MM/PBSA-style calculations
- `solvateoct` for a truncated octahedral water box

### Input: `addWater.in`

```bash
#!/bin/bash

source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.opc

loadamberprep   file.prepin
loadamberparams file.frcmod

set default PBradii mbondi3

w = loadpdb minH.pdb

# --- Disulfide bonds ---
bond w.7.SG w.33.SG
bond w.19.SG w.130.SG
bond w.77.SG w.87.SG

addIons w Na+ 1
solvateoct w OPCBOX 12.0

check w
charge w

savepdb w addWAT.pdb
saveamberparm w addWAT.top addWAT.crd

quit
```

### Run

```bash
tleap -f addWater.in
```

### Outputs

```text
addWAT.pdb
addWAT.top
addWAT.crd
```

### Important Ion Note

`addIons w Na+ 1` adds exactly one sodium ion.

If your goal is neutralization, consider checking the system charge and using:

```bash
addIons w Na+ 0
addIons w Cl- 0
```

Use `charge w` inside `tleap` to verify the final charge.

---

## 6. Stage 4 — `minWAT`: Minimize Water/Ions with Solute Restrained

### Purpose

`minWAT` relaxes water and ions while keeping the enzyme/solute fixed or strongly restrained.

### Input: `minWAT.in`

```text
Minimizatrion of water molecules ---Fixing Enzyme---
 &cntrl
 imin   = 1,
 maxcyc = 3000,
 ncyc   = 1500,
 ntb    = 1,
 ntr    = 1,
 cut    = 12
 /
Hold the ENZYME fixed
500.0
RES 1 302
END
END
```

### Run: `minWAT.job`

```bash
#!/bin/bash

mpirun -np 10 pmemd.MPI -O \
  -i minWAT.in \
  -o minWAT.out \
  -p addWAT.top \
  -c addWAT.crd \
  -r minWAT.restrt \
  -ref addWAT.crd
```

### Output

```text
minWAT.restrt
minWAT.out
```

---

## 7. Stage 5 — `minALL`: Minimize Full Solvated System

### Purpose

`minALL` performs full-system minimization after water relaxation.

### Input: `minALL.in`

```text
Minimizatrion of water molecules ---Fixing Enzyme---
 &cntrl
 imin   = 1,
 maxcyc = 3000,
 ncyc   = 1500,
 ntb    = 1,
 ntr    = 0,
 cut    = 12
 /
Hold the ENZYME fixed
500.0
END
END
```

### Cleaner Recommended Version

Because `ntr=0`, the restraint block is not needed. A cleaner file is:

```text
Full-system minimization
 &cntrl
 imin   = 1,
 maxcyc = 3000,
 ncyc   = 1500,
 ntb    = 1,
 ntr    = 0,
 cut    = 12.0,
 /
```

### Run: `minALL.job`

```bash
#!/bin/bash

mpirun -np 10 $PMEMDHOME/bin/pmemd.MPI -O \
  -i minALL.in \
  -o minALL.out \
  -p addWAT.top \
  -c minWAT.restrt \
  -r min_all.restrt \
  -ref minWAT.restrt
```

### Output

```text
min_all.restrt
minALL.out
```

---

## 8. Stage 6 — Staged Minimization, Heating, and Equilibration

After `min_all.restrt`, your workflow runs:

```text
min1 → min2 → min3 → min4 → heat → eq1 → eq2
```

using `min_heat.qsub`.

### Input: `min_heat.qsub`

```bash
#!/bin/bash

cp min_all.restrt complex_ini.rst7

old=ini

for name in min1 min2 min3 min4 heat eq1 eq2 ; do

 $PMEMDHOME/bin/pmemd.cuda -O -i $name.i -o complex_$name.mdout -p addWAT.top -c complex_$old.rst7 \
 -ref min_all.restrt -x complex_$name.nc -r complex_$name.rst7 -inf complex_$name.mdinfo
  old=$name

done
```

### How the Restart Chain Works

The script first makes:

```bash
cp min_all.restrt complex_ini.rst7
```

Then the loop begins:

```text
old=ini
name=min1
```

This means AMBER reads:

```text
complex_ini.rst7
```

and writes:

```text
complex_min1.rst7
```

Next:

```text
old=min1
name=min2
```

It reads:

```text
complex_min1.rst7
```

and writes:

```text
complex_min2.rst7
```

This continues until:

```text
complex_eq2.rst7
```

The final equilibration file becomes the starting point for production MD.

---

## 9. Stage 6A — `min1.i`: Minimize Water, Hydrogens, and Ions with Solute Restraints

### Purpose

`min1.i` relaxes water, hydrogen atoms, and ions while applying strong restraints to the solute.

### Input: `min1.i`

```text
Minimization of water, hydrogens atoms and ions with restraints for solute
&cntrl
 imin=1,                         ! Turn on minimization
 ncyc=1000,                      ! Number of steepest descent steps
 maxcyc=10000,                   ! Total number of minimization cycles
 ntmin=1,                        ! Steepest descent for ncyc steps, then conjugate gradient
 ntx=1,                          ! Coordinates, but no velocities, will be read
 cut=12.0,                       ! 12 A cutoff required for OPC water
 ntwx=500,                       ! Coordinates written every ntwx steps
 ntpr=50,                        ! Print out energy information every ntpr steps
 ioutfm=1,                       ! Use Binary NetCDF trajectory format
 iwrap=0,                        ! No wrapping
 ntxo=1,                         ! NetCDF file
 ntr=1,                          ! Restraints on
 restraint_wt=25.0,              ! kcal/mol/A**2 restraint force constant
 restraintmask='!(:WAT,Na+,Cl-  | @H=)',
 nmropt=0
/
```

### Key Parameters

| Parameter | Meaning |
|---|---|
| `imin=1` | Run minimization |
| `maxcyc=10000` | Total minimization cycles |
| `ncyc=1000` | Steepest descent before conjugate gradient |
| `ntb=1` | Periodic boundary conditions at constant volume |
| `cut=12.0` | Nonbonded cutoff, consistent with OPC water |
| `ntr=1` | Positional restraints active |
| `restraint_wt=25.0` | Strong restraint |
| `restraintmask='!(:WAT,Na+,Cl-  | @H=)'` | Restrain all atoms except water, ions, and hydrogens |

### Interpretation

This stage allows the most mobile and newly placed atoms to relax first while keeping the heavier solute framework stable.

---

## 10. Stage 6B — `min2.i`: Minimize Side Chains

### Purpose

`min2.i` continues minimization while allowing side-chain relaxation and keeping backbone-like atoms restrained.

### Input: `min2.i`

```text
Minimization of side chains
&cntrl
 imin=1,                         ! Turn on minimization
 ncyc=1000,                      ! Number of steepest descent steps
 maxcyc=10000,                   ! Total number of minimization cycles
 ntmin=1,                        ! Steepest descent for ncyc steps, then conjugate gradient
 ntx=1,                          ! Coordinates, but no velocities, will be read
 cut=12.0,                       ! 12 A cutoff required for OPC water
 ntwx=500,                       ! Coordinates written every ntwx steps
 ntpr=50,                        ! Print out energy information every ntpr steps
 ioutfm=1,                       ! Use Binary NetCDF trajectory format
 iwrap=0,                        ! No wrapping
 ntxo=1,                         ! NetCDF file
 ntr=1,                          ! Restraints on
 restraint_wt=25.0,              ! kcal/mol/A**2 restraint force constant
 restraintmask='!:WAT&@CA,C,O,N,ZN',
 nmropt=0
/
```

### Key Restraint Mask

```text
restraintmask='!:WAT&@CA,C,O,N,ZN'
```

This mask is intended to restrain selected non-water solute atoms such as backbone atoms and `ZN`.

> **Check carefully:** If your system does not contain `ZN`, the mask may still work, but it is cleaner to remove `ZN`. If your system contains another metal or no metal, update the mask.

---

## 11. Stage 6C — `min3.i`: Minimize Everything Except CA

### Purpose

`min3.i` relaxes most of the system while keeping alpha carbons and zinc restrained.

### Input: `min3.i`

```text
Minimization of everything but CA
&cntrl
 imin=1,                         ! Turn on minimization
 ncyc=1000,                      ! Number of steepest descent steps
 maxcyc=10000,                   ! Total number of minimization cycles
 ntmin=1,                        ! Steepest descent for ncyc steps, then conjugate gradient
 ntx=1,                          ! Coordinates, but no velocities, will be read
 cut=12.0,                       ! 12 A cutoff required for OPC water
 ntwx=500,                       ! Coordinates written every ntwx steps
 ntpr=50,                        ! Print out energy information every ntpr steps
 ioutfm=1,                       ! Use Binary NetCDF trajectory format
 iwrap=0,                        ! No wrapping
 ntxo=1,                         ! NetCDF file
 ntr=1,                          ! Restraints on
 restraint_wt=25.0,              ! kcal/mol/A**2 restraint force constant
 restraintmask='@CA,ZN',
 nmropt=0
/
```

### Key Restraint Mask

```text
restraintmask='@CA,ZN'
```

This keeps `CA` atoms and `ZN` atoms restrained while allowing other atoms to relax.

> **Check carefully:** If your protein does not contain a zinc ion, remove `ZN` from the mask.

---

## 12. Stage 6D — `min4.i`: Full-System Minimization Without Restraints

### Purpose

`min4.i` performs final minimization without restraints.

### Input: `min4.i`

```text
Minimization of everything
&cntrl
 imin=1,                         ! Turn on minimization
 ncyc=1000,                      ! Number of steepest descent steps
 maxcyc=10000,                   ! Total number of minimization cycles
 ntmin=1,                        ! Steepest descent for ncyc steps, then conjugate gradient
 ntx=1,                          ! Coordinates, but no velocities, will be read
 cut=12.0,                       ! 12 A cutoff required for OPC water
 ntwx=500,                       ! Coordinates written every ntwx steps
 ntpr=50,                        ! Print out energy information every ntpr steps
 ioutfm=1,                       ! Use Binary NetCDF trajectory format
 iwrap=0,                        ! No wrapping
 ntxo=1,                         ! NetCDF file
 ntr=0,                          ! No restraints
 nmropt=0
/
```

### Key Parameter

```text
ntr=0
```

This means no positional restraints are active.

---

## 13. Stage 6E — `heat.i`: Heating from 10 K to 310 K

### Purpose

`heat.i` gradually heats the system from 10 K to 310 K over 200 ps under NVT conditions with backbone restraints.

### Input: `heat.i`

```text
MD heating of system over 200 ps
&cntrl
 imin=0,                         ! Not a minimisation run
 irest=0,                        ! New simulation
 ntx=1,                          ! Read coordinates but not velocities
 nscm=1000,                      ! Reset COM every 1000 steps
 nstlim=100000, dt=0.002,        ! Run MD for 200 ps with a timestep of 2 fs
 ntpr=500, ntwx=5000,            ! Write trajectory every 10 ps, energies every 1 ps
 ioutfm=1,                       ! Use Binary NetCDF trajectory format
 iwrap=0,                        ! No wrapping
 ntxo=1,                         ! NetCDF file
 ntwr=5000,                      ! Write restart file every 10 ps
 ntb=1,                          ! PBC at constant volume (NVT)
 cut=12.0,                       ! 12 A cutoff required for OPC water
 ntc=2, ntf=2,                   ! SHAKE on all H
 ntp=0,                          ! No pressure regulation
 ntt=3,                          ! Langevin thermostat
 gamma_ln=5.0,                   ! Langevin collision frequency
 ig=-1,                          ! Randomize RNG seed
 ntr=1,                          ! Positional restraints on
 restraint_wt=5.0,               ! kcal/mol/A**2 restraint force constant
 restraintmask='!:WAT&@CA,N,C,O,ZN',
 nmropt=1,                       ! NMR restraints for temperature ramp
/
&wt
  TYPE='TEMP0', ISTEP1=0, ISTEP2=100000,
  VALUE1=10.0, VALUE2=310.0,
/
&wt TYPE='END' /
```

### Simulation Time

```text
nstlim × dt = 100000 × 0.002 ps = 200 ps
```

### Temperature Ramp

The temperature is controlled using `nmropt=1` and a `&wt` block:

```text
&wt
  TYPE='TEMP0', ISTEP1=0, ISTEP2=100000,
  VALUE1=10.0, VALUE2=310.0,
/
&wt TYPE='END' /
```

This means the target temperature increases from:

```text
10 K → 310 K
```

over 100,000 MD steps.

### Key Parameters

| Parameter | Meaning |
|---|---|
| `imin=0` | Run MD, not minimization |
| `irest=0` | New MD simulation |
| `ntx=1` | Read coordinates only, no velocities |
| `ntb=1` | Constant volume periodic simulation |
| `ntp=0` | No pressure coupling |
| `ntt=3` | Langevin thermostat |
| `gamma_ln=5.0` | Collision frequency |
| `ntc=2, ntf=2` | SHAKE on bonds involving hydrogen |
| `ntr=1` | Positional restraints active |
| `restraint_wt=5.0` | Backbone restraint strength |
| `restraintmask='!:WAT&@CA,N,C,O,ZN'` | Restrain selected solute atoms |
| `ig=-1` | Randomize seed |

### Heating Quality Checks

After heating, check:

```bash
tail -n 50 complex_heat.mdout
```

Look for:

- temperature gradually reaches 310 K
- no SHAKE failure
- no abnormal energy jump
- no `NaN`
- structure remains stable

---

## 14. Stage 6F — `eq1.i`: 500 ps NPT Equilibration with Backbone Restraints

### Purpose

`eq1.i` equilibrates the system under NPT conditions while maintaining backbone restraints.

### Input: `eq1.i`

```text
500 ps NPT MD equilibration with backbone restraints
&cntrl
 imin=0,                         ! Not a minimisation run
 irest=1,                        ! Restart — carry velocities from heat
 ntx=5,                          ! Read coordinates AND velocities
 nscm=1000,                      ! Reset COM every 1000 steps
 nstlim=250000, dt=0.002,        ! Run MD for 500 ps with a timestep of 2 fs
 ntpr=500, ntwx=5000,            ! Write trajectory every 10 ps, energies every 1 ps
 ioutfm=1,                       ! Use Binary NetCDF trajectory format
 iwrap=0,                        ! No wrapping
 ntxo=1,                         ! NetCDF file
 ntwr=250000,                    ! Write restart at end of run
 ntb=2,                          ! PBC at constant pressure (NPT)
 cut=12.0,                       ! 12 A cutoff required for OPC water
 ntc=2, ntf=2,                   ! SHAKE on all H
 ntp=1,                          ! Isotropic pressure regulation
 pres0=1.01325,                  ! Reference pressure in bars
 taup=2.0,                       ! Softer pressure relaxation time (2 ps — safer for GPU NPT)
 barostat=2,                     ! Monte Carlo barostat — more stable than Berendsen for GPU
 ntt=3,                          ! Langevin thermostat
 gamma_ln=2.0,                   ! Langevin collision frequency
 tempi=310.0,                    ! Initial temperature in K
 temp0=310.0,                    ! Target temperature in K
 ig=-1,                          ! Randomize RNG seed
 ntr=1,                          ! Positional restraints on
 restraint_wt=5.0,               ! kcal/mol/A**2 restraint force constant
 restraintmask='!:WAT&@CA,N,C,O,ZN',
 nmropt=0,
/
```

### Simulation Time

```text
nstlim × dt = 250000 × 0.002 ps = 500 ps
```

### Key Parameters

| Parameter | Meaning |
|---|---|
| `irest=1, ntx=5` | Continue from heating with velocities |
| `ntb=2` | Constant pressure periodic simulation |
| `ntp=1` | Isotropic pressure regulation |
| `barostat=2` | Monte Carlo barostat |
| `pres0=1.01325` | Reference pressure in bar |
| `taup=2.0` | Pressure relaxation time |
| `temp0=310.0` | Target temperature |
| `ntr=1` | Positional restraints active |
| `restraint_wt=5.0` | Strong backbone restraint |

### Equilibration Checks

Check:

```bash
tail -n 50 complex_eq1.mdout
```

Look for:

- stable temperature around 310 K
- pressure fluctuations are normal
- density begins to stabilize
- volume adjusts reasonably
- no SHAKE or `vlimit` errors

---

## 15. Stage 6G — `eq2.i`: 500 ps NPT Equilibration with Reduced Restraints

### Purpose

`eq2.i` continues NPT equilibration with weaker restraints.

### Input: `eq2.i`

```text
500 ps NPT MD equilibration with reduced backbone restraints
&cntrl
 imin=0,                         ! Not a minimisation run
 irest=1,                        ! Restart simulation
 ntx=5,                          ! Read coordinates AND velocities
 nscm=1000,                      ! Reset COM every 1000 steps
 nstlim=250000, dt=0.002,        ! Run MD for 500 ps with a timestep of 2 fs
 ntpr=500, ntwx=5000,            ! Write trajectory every 10 ps, energies every 1 ps
 ioutfm=1,                       ! Use Binary NetCDF trajectory format
 iwrap=0,                        ! No wrapping
 ntxo=1,                         ! NetCDF file
 ntwr=250000,                    ! Write restart at end of run
 ntb=2,                          ! PBC at constant pressure (NPT)
 cut=12.0,                       ! 12 A cutoff required for OPC water
 ntc=2, ntf=2,                   ! SHAKE on all H
 ntp=1,                          ! Isotropic pressure regulation
 pres0=1.01325,                  ! Reference pressure in bars
 taup=2.0,                       ! Softer pressure relaxation time (2 ps)
 barostat=2,                     ! Monte Carlo barostat — more stable than Berendsen for GPU
 ntt=3,                          ! Langevin thermostat
 gamma_ln=2.0,                   ! Langevin collision frequency
 tempi=310.0,                    ! Initial temperature in K
 temp0=310.0,                    ! Target temperature in K
 ig=-1,                          ! Randomize RNG seed
 ntr=1,                          ! Positional restraints on
 restraint_wt=2.0,               ! Reduced to 2.0 kcal/mol/A**2 (was 5.0 in eq1)
 restraintmask='!:WAT&@CA,N,C,O,ZN',
 nmropt=0,
/
```

### Simulation Time

```text
250000 × 0.002 ps = 500 ps
```

### Key Difference from `eq1.i`

```text
restraint_wt=2.0
```

This reduces the restraint strength from 5.0 to 2.0 kcal/mol/Å², allowing more relaxation while still protecting the backbone from abrupt movement.

### Final Output

After the loop, the most important file is:

```text
complex_eq2.rst7
```

This file should be used as the starting coordinate/restart file for production MD.

---

## 16. Run the Full Staged Min/Heat/Eq Workflow

Run:

```bash
bash min_heat.qsub
```

or, if submitting to a cluster scheduler, adapt the command based on your system.

After completion, check:

```bash
ls complex_min1.rst7
ls complex_min2.rst7
ls complex_min3.rst7
ls complex_min4.rst7
ls complex_heat.rst7
ls complex_eq1.rst7
ls complex_eq2.rst7
```

Also check output files:

```bash
tail -n 30 complex_min1.mdout
tail -n 30 complex_min2.mdout
tail -n 30 complex_min3.mdout
tail -n 30 complex_min4.mdout
tail -n 30 complex_heat.mdout
tail -n 30 complex_eq1.mdout
tail -n 30 complex_eq2.mdout
```

---

## 17. Production MD After Equilibration

After `complex_eq2.rst7` is generated from the equilibration stage, production MD can begin.

Your production step uses two main files:

```text
runMD.qsub
md.i
```

The production workflow starts from:

```text
complex_eq2.rst7
```

and generates sequential production chunks:

```text
complex_md1.nc,  complex_md1.rst7
complex_md2.nc,  complex_md2.rst7
complex_md3.nc,  complex_md3.rst7
...
complex_md50.nc, complex_md50.rst7
```

This chunked-restart strategy is useful because long simulations can be continued safely from the last completed restart file.

---

## 17.1 Production Input File: `md.i`

Your `md.i` file defines each production MD chunk as **2 ns NPT production MD without restraints**.

```text
2 ns NPT production MD — no restraints
&cntrl
 imin=0,                         ! Not a minimisation run
 irest=1,                        ! Restart simulation
 ntx=5,                          ! Read coordinates AND velocities
 nscm=1000,                      ! Reset COM every 1000 steps
 nstlim=1000000, dt=0.002,       ! Run MD for 2 ns with a timestep of 2 fs
 ntpr=500, ntwx=5000,            ! Write trajectory every 10 ps, energies every 1 ps
 ioutfm=1,                       ! Use Binary NetCDF trajectory format
 iwrap=1,                        ! Wrap coordinates into primary box
 ntxo=1,                         ! NetCDF file
 ntwr=250000,                    ! Write restart every 500 ps
 ntb=2,                          ! PBC at constant pressure (NPT)
 cut=12.0,                       ! 12 A cutoff required for OPC water
 ntc=2, ntf=2,                   ! SHAKE on all H
 ntp=1,                          ! Isotropic pressure regulation
 pres0=1.01325,                  ! Reference pressure in bars
 taup=2.0,                       ! Softer pressure relaxation time (2 ps)
 barostat=2,                     ! Monte Carlo barostat — more stable than Berendsen for GPU
 ntt=3,                          ! Langevin thermostat
 gamma_ln=2.0,                   ! Langevin collision frequency
 tempi=310.0,                    ! Initial temperature in K
 temp0=310.0,                    ! Target temperature in K
 ig=-1,                          ! Randomize RNG seed
 ntr=0,                          ! No positional restraints in production
 nmropt=0,
/
```

### Key Parameters in `md.i`

| Parameter | Meaning |
|---|---|
| `imin=0` | Run molecular dynamics, not minimization |
| `irest=1` | Restart from a previous simulation |
| `ntx=5` | Read coordinates and velocities from the restart file |
| `nstlim=1000000` | Run 1,000,000 MD steps |
| `dt=0.002` | 2 fs timestep |
| `ntpr=500` | Print energy information every 500 steps |
| `ntwx=5000` | Write trajectory every 5000 steps |
| `ntwr=250000` | Write restart every 250,000 steps |
| `ntb=2` | Periodic boundary conditions with constant pressure |
| `ntp=1` | Isotropic pressure coupling |
| `barostat=2` | Monte Carlo barostat |
| `ntt=3` | Langevin thermostat |
| `gamma_ln=2.0` | Langevin collision frequency |
| `temp0=310.0` | Target temperature is 310 K |
| `ntc=2, ntf=2` | SHAKE constraints for bonds involving hydrogen |
| `ntr=0` | No positional restraints during production |
| `ioutfm=1` | Write NetCDF trajectory |
| `iwrap=1` | Wrap coordinates into the primary box |

### Simulation Time Per Chunk

```text
nstlim × dt = 1,000,000 × 0.002 ps = 2,000 ps = 2 ns
```

Therefore, each production file such as `complex_md1.nc` represents approximately **2 ns** of simulation.

---

## 17.2 Production Running Script: `runMD.qsub`

Your active production script starts from `eq2` and runs `md1` to `md50`.

```bash
#!/bin/bash

old=eq2
for name in md1 md2 md3 md4 md5 md6 md7 md8 md9 md10 md11 md12 md13 md14 md15 md16 md17 md18 md19 md20 md21 md22 md23 md24 \
md25 md26 md27 md28 md29 md30 md31 md32 md33 md34 md35 md36 md37 md38 md39 md40 md41 md42 md43 md44 md45 md46 md47 md48 md49 md50 ; do

#old=md50
#for name in md51 md52 md53 md54 md55 md56 md57 md58 md59 md60 md61 md62 md63 md64 md65 md66 md67 md68 md69 md70 \
#md71 md72 md73 md74 md75 md76 md77 md78 md79 md80 md81 md82 md83 md84 md85 md86 md87 md88 md89 md90 md91 md92 md93 md94 md95 md96 md97 md98 md99 md100 \
#md101 md102 md103 md104 md105 md106 md107 md108 md109 md110 md111 md112 md113 md114 md115 md116 md117 md118 md119 md120 md121 md122 md123 md124 md125 \
#md126 md127 md128 md129 md130 md131 md132 md133 md134 md135 md136 md137 md138 md139 md140 md141 md142 md143 md144 md145 md146 md147 md148 md149 md150 \
#md151 md152 md153 md154 md155 md156 md157 md158 md159 md160 md161 md162 md163 md164 md165 md166 md167 md168 md169 md170 md171 md172 md173 md174 md175 \
#md176 md177 md178 md179 md180 md181 md182 md183 md184 md185 md186 md187 md188 md189 md190 md191 md192 md193 md194 md195 md196 md197 md198 md199 md200 \
#md201 md202 md203 md204 md205 md206 md207 md208 md209 md210 md211 md212 md213 md214 md215 md216 md217 md218 md219 md220 md221 md222 md223 md224 md225 \
#md226 md227 md228 md229 md230 md231 md232 md233 md234 md235 md236 md237 md238 md239 md240 md241 md242 md243 md244 md245 md246 md247 md248 md249 md250 ; do

#old=md250
#for name in md251 md252 md253 md254 md255 md256 md257 md258 md259 md260 md261 md262 md263 md264 md265 md266 md267 md268 md269 md270 \
#md271 md272 md273 md274 md275 md276 md277 md278 md279 md280 md281 md282 md283 md284 md285 md286 md287 md288 md289 md290 \
#md291 md292 md293 md294 md295 md296 md297 md298 md299 md300 md301 md302 md303 md304 md305 md306 md307 md308 md309 md310 \
#md311 md312 md313 md314 md315 md316 md317 md318 md319 md320 md321 md322 md323 md324 md325 md326 md327 md328 md329 md330 \
#md331 md332 md333 md334 md335 md336 md337 md338 md339 md340 md341 md342 md343 md344 md345 md346 md347 md348 md349 md350 \
#md351 md352 md353 md354 md355 md356 md357 md358 md359 md360 md361 md362 md363 md364 md365 md366 md367 md368 md369 md370 \
#md371 md372 md373 md374 md375 md376 md377 md378 md379 md380 md381 md382 md383 md384 md385 md386 md387 md388 md389 md390 \
#md391 md392 md393 md394 md395 md396 md397 md398 md399 md400 md401 md402 md403 md404 md405 md406 md407 md408 md409 md410 \
#md411 md412 md413 md414 md415 md416 md417 md418 md419 md420 md421 md422 md423 md424 md425 md426 md427 md428 md429 md430 \
#md431 md432 md433 md434 md435 md436 md437 md438 md439 md440 md441 md442 md443 md444 md445 md446 md447 md448 md449 md450 \
#md451 md452 md453 md454 md455 md456 md457 md458 md459 md460 md461 md462 md463 md464 md465 md466 md467 md468 md469 md470 \
#md471 md472 md473 md474 md475 md476 md477 md478 md479 md480 md481 md482 md483 md484 md485 md486 md487 md488 md489 md490 \
#md491 md492 md493 md494 md495 md496 md497 md498 md499 md500 ; do

 $PMEMDHOME/bin/pmemd.cuda -O -i md.i -o complex_$name.mdout -p addWAT.top -c complex_$old.rst7 \
 -ref min_all.restrt -x complex_$name.nc -r complex_$name.rst7 -inf complex_$name.mdinfo
 old=$name
 done
```

### Active Production Block

The active part is:

```bash
old=eq2

for name in md1 md2 md3 md4 md5 md6 md7 md8 md9 md10 md11 md12 md13 md14 md15 md16 md17 md18 md19 md20 md21 md22 md23 md24 \
md25 md26 md27 md28 md29 md30 md31 md32 md33 md34 md35 md36 md37 md38 md39 md40 md41 md42 md43 md44 md45 md46 md47 md48 md49 md50 ; do

 $PMEMDHOME/bin/pmemd.cuda -O -i md.i -o complex_$name.mdout -p addWAT.top -c complex_$old.rst7 \
 -ref min_all.restrt -x complex_$name.nc -r complex_$name.rst7 -inf complex_$name.mdinfo

 old=$name
done
```

### How the Restart Chain Works

At the first production step:

```text
old=eq2
name=md1
```

AMBER reads:

```text
complex_eq2.rst7
```

and writes:

```text
complex_md1.rst7
complex_md1.nc
complex_md1.mdout
complex_md1.mdinfo
```

Then the script sets:

```bash
old=md1
```

At the second production step:

```text
old=md1
name=md2
```

AMBER reads:

```text
complex_md1.rst7
```

and writes:

```text
complex_md2.rst7
complex_md2.nc
```

This continues until:

```text
complex_md50.rst7
complex_md50.nc
```

---

## 17.3 Total Simulation Time

Because each chunk is 2 ns:

| Production Range | Number of Chunks | Total Time |
|---|---:|---:|
| `md1`–`md50` | 50 | 100 ns |
| `md51`–`md250` | 200 | 400 ns additional |
| `md1`–`md250` | 250 | 500 ns total |
| `md251`–`md500` | 250 | 500 ns additional |
| `md1`–`md500` | 500 | 1000 ns = 1 µs total |

Your script already contains commented continuation blocks for:

```text
md51 → md250
md251 → md500
```

To continue the simulation, uncomment the correct continuation block and set `old` to the last completed restart file.

Example:

```bash
old=md50
for name in md51 md52 md53 ... md250 ; do
   ...
done
```

---

## 17.4 Cleaner Optional Production Script Using `seq`

Your original script is valid, but the long manual list is easy to mistype. A cleaner optional version is:

```bash
#!/bin/bash

start=1
end=50
old=eq2

for i in $(seq $start $end); do
    name=md${i}

    $PMEMDHOME/bin/pmemd.cuda -O \
      -i md.i \
      -o complex_${name}.mdout \
      -p addWAT.top \
      -c complex_${old}.rst7 \
      -ref min_all.restrt \
      -x complex_${name}.nc \
      -r complex_${name}.rst7 \
      -inf complex_${name}.mdinfo

    old=${name}
done
```

To continue from `md50` to `md250`, use:

```bash
start=51
end=250
old=md50
```

To continue from `md250` to `md500`, use:

```bash
start=251
end=500
old=md250
```

---

## 17.5 Files Required Before Production

Before running `runMD.qsub`, check:

```bash
ls addWAT.top
ls complex_eq2.rst7
ls min_all.restrt
ls md.i
ls runMD.qsub
```

The most important file is:

```text
complex_eq2.rst7
```

If this file does not exist, production cannot start because the first production chunk uses:

```bash
-c complex_eq2.rst7
```

---

## 17.6 Run Production MD

Run locally:

```bash
bash runMD.qsub
```

If your cluster uses a scheduler, submit according to your HPC system, for example:

```bash
qsub runMD.qsub
```

or adapt it for SLURM:

```bash
sbatch runMD.sh
```

depending on the cluster.

---

## 17.7 Production Output Files

Each production chunk produces:

| Output File | Meaning |
|---|---|
| `complex_mdX.mdout` | AMBER production output log |
| `complex_mdX.nc` | NetCDF trajectory for production chunk X |
| `complex_mdX.rst7` | Restart file for the next chunk |
| `complex_mdX.mdinfo` | Runtime status information |

Example for `md1`:

```text
complex_md1.mdout
complex_md1.nc
complex_md1.rst7
complex_md1.mdinfo
```

---

## 17.8 Production Monitoring Commands

### Check if files are being produced

```bash
ls complex_md1.*
```

### Monitor runtime status

```bash
tail -f complex_md1.mdinfo
```

### Check final output of a chunk

```bash
tail -n 50 complex_md1.mdout
```

### Check the last completed restart

```bash
ls complex_md*.rst7
```

### Count trajectory chunks

```bash
ls complex_md*.nc | wc -l
```

---

## 17.9 Production Quality-Control Checks

After each chunk or after a group of chunks, check:

```bash
tail -n 30 complex_md1.mdout
tail -n 30 complex_md50.mdout
```

Look for:

- normal termination
- no `NaN`
- no `SHAKE failure`
- no `vlimit exceeded`
- temperature close to 310 K
- reasonable pressure fluctuations
- no sudden energy explosion
- restart file was written successfully

### Important Values to Watch

| Quantity | What to Check |
|---|---|
| Temperature | Should fluctuate around `temp0=310.0` K |
| Pressure | Can fluctuate strongly in NPT but should not diverge |
| Density | Should be physically reasonable and stable after equilibration |
| Total energy | Should not show sudden nonphysical jumps |
| Restart files | Must exist for continuation |
| Trajectory files | Must be readable by `cpptraj` |

---

## 17.10 Common Production Errors

### Error 1: `complex_eq2.rst7` Not Found

Cause:

The equilibration step did not finish or the file name is different.

Fix:

```bash
ls complex_eq2.rst7
```

If missing, rerun or check the `eq2` stage.

---

### Error 2: `complex_md50.rst7` Missing When Continuing to `md51`

Cause:

The `md50` chunk did not finish, or the file is in another folder.

Fix:

```bash
ls complex_md50.rst7
tail -n 50 complex_md50.mdout
```

Continue from the last completed restart file.

---

### Error 3: SHAKE Failure During Production

Likely causes:

- insufficient equilibration
- unstable ligand parameters
- bad contacts
- too aggressive timestep for the system
- unstable water/ion placement
- metal or nonstandard residue parameter problem

Possible fixes:

```text
1. Inspect the final equilibration structure.
2. Rerun minimization/equilibration with stronger restraints.
3. Temporarily reduce dt from 0.002 to 0.001.
4. Check ligand parameters and charge.
5. Check for unusual atoms, metals, or covalent groups.
```

---

### Error 4: Production Runs but Trajectory Looks Broken

Likely cause:

- Periodic boundary imaging artifacts
- Molecules wrapped across the box

Fix:

Use `cpptraj` with:

```text
autoimage
center
image
```

Do not interpret raw wrapped trajectories without processing.

---

## 17.11 Method Description for Production MD

```text
Production molecular dynamics was performed using the GPU-accelerated AMBER engine pmemd.cuda. The final equilibrated restart file, complex_eq2.rst7, was used as the starting structure. Production simulations were carried out under NPT conditions at 310 K and 1.01325 bar using a 2 fs timestep, Langevin temperature control, Monte Carlo barostat, SHAKE constraints on bonds involving hydrogen, and a 12 Å nonbonded cutoff. Each production segment was run for 2 ns, producing NetCDF trajectory files and restart files for continuation. Sequential restart chaining was used, where each production segment used the restart file generated from the previous segment.
```

---

## 17.12 Final Production Rule

Never start production MD until these files are confirmed:

```bash
ls addWAT.top
ls min_all.restrt
ls complex_eq2.rst7
ls md.i
```

The production chain depends on this sequence:

```text
complex_eq2.rst7
   ↓ md1
complex_md1.rst7
   ↓ md2
complex_md2.rst7
   ↓ md3
...
```

If one restart file is missing, all later production chunks will fail.




## 18. Quality-Control Checklist for Heat and Equilibration

### After `min1`–`min4`

Check:

- energy decreases
- no `NaN`
- no abnormal geometry
- no minimization failure
- restraints are being applied as expected

### After Heating

Check:

- temperature reaches 310 K gradually
- total energy does not explode
- no SHAKE failure
- no `vlimit exceeded`
- no sudden structural distortion

### After `eq1`

Check:

- temperature stable around 310 K
- pressure fluctuates but does not diverge
- density starts stabilizing
- no water box collapse
- no ligand ejection unless scientifically expected

### After `eq2`

Check:

- system remains stable with weaker restraints
- density and volume are reasonable
- structure is suitable for production
- `complex_eq2.rst7` exists

---

## 19. Common Problems and Fixes

### Problem 1: `complex_ini.rst7` Missing

Cause:

`min_heat.qsub` expects to create this file from:

```bash
cp min_all.restrt complex_ini.rst7
```

If `min_all.restrt` does not exist, the workflow fails.

Fix:

```bash
ls min_all.restrt
```

If missing, rerun `minALL.job`.

---

### Problem 2: `pmemd.cuda` Cannot Find Input File

Cause:

The loop expects files named:

```text
min1.i
min2.i
min3.i
min4.i
heat.i
eq1.i
eq2.i
```

Fix:

If your uploaded files are named like `min1(1).i`, rename them:

```bash
mv 'min1(1).i' min1.i
mv 'min2(1).i' min2.i
mv 'min3(1).i' min3.i
mv 'min4(1).i' min4.i
mv 'heat(1).i' heat.i
mv 'eq1(1).i' eq1.i
mv 'eq2(1).i' eq2.i
mv 'min_heat(1).qsub' min_heat.qsub
```

---

### Problem 3: Restraint Mask Contains `ZN` but System Has No Zinc

Cause:

Several masks include:

```text
ZN
```

If your system has no zinc, it may be unnecessary or confusing.

Fix:

Check:

```bash
grep ZN addWAT.pdb
```

If no zinc is present, remove `ZN` from masks.

---

### Problem 4: Heating Fails with SHAKE Error

Likely causes:

- minimization was insufficient
- bad contacts remain
- ligand parameters are poor
- timestep is too large for unstable starting geometry
- heating is too fast

Fixes:

- inspect `complex_min4.rst7`
- run additional minimization
- temporarily reduce `dt=0.001`
- reduce heating speed
- check ligand `.frcmod`
- verify no overlapping atoms

---

### Problem 5: NPT Equilibration Fails

Likely causes:

- unstable system after heating
- pressure coupling begins too early
- bad density or box size
- poor solvation geometry
- bad ligand or ion parameters

Fixes:

- extend NVT heating/equilibration
- use stronger restraints initially
- inspect density/volume in `mdout`
- visualize the structure
- repeat water minimization if needed

---

## 20. Complete Command Order

Use this as the complete clean command order:

```bash
# Stage 1: dry complex
tleap -f addH.in

# Stage 2: hydrogen minimization
bash minH.job

# Required bridge step
ambpdb -p complex_addH.top -c complex_minH.restrt > minH.pdb

# Stage 3: add water and ions
tleap -f addWater.in

# Stage 4: solvent minimization
bash minWAT.job

# Stage 5: full minimization
bash minALL.job

# Stage 6: staged minimization, heating, and equilibration
bash min_heat.qsub

# Stage 7: production MD
# Option 1: run locally
bash runMD.qsub

# Option 2: submit to a PBS/Torque cluster, if applicable
qsub runMD.qsub

# Stage 8: trajectory analysis
bash cpptraj-analysis.sh
```

---

## 21. Method Description Template

```text
The protein-ligand complex was first processed using tleap with the ff19SB protein force field and GAFF2 ligand parameters. Ligand topology and parameter files generated from antechamber and parmchk2 were loaded as file.prepin and file.frcmod. Disulfide bonds were manually assigned, and the dry hydrogen-added complex was saved as complex_addH.top and complex_addH.crd. Hydrogen atoms were minimized using sander, and the minimized restart was converted to minH.pdb using ambpdb.

The hydrogen-minimized complex was then solvated in an OPC truncated octahedral water box using tleap, with ions added as required. The solvated system was saved as addWAT.top and addWAT.crd. Water and ions were minimized while the solute was restrained, followed by full-system minimization.

Additional staged minimizations were performed using decreasing restraint strategies, followed by gradual heating from 10 K to 310 K over 200 ps under NVT conditions. The system was then equilibrated under NPT conditions for 500 ps with backbone restraints and another 500 ps with reduced restraints. The final equilibrated restart file, complex_eq2.rst7, was used as the starting structure for production molecular dynamics.
```

---

## Final Recommendation

Your current heat/equilibration stage is logically structured:

```text
min1 → min2 → min3 → min4 → heat → eq1 → eq2
```

The most important practical checks are:

1. Make sure the uploaded filenames are renamed to match the loop names: `min1.i`, `min2.i`, `min3.i`, `min4.i`, `heat.i`, `eq1.i`, and `eq2.i`.
2. Make sure `min_all.restrt` exists before running `min_heat.qsub`.
3. Confirm whether your system really contains `ZN`; if not, remove `ZN` from restraint masks.
4. Confirm that `complex_eq2.rst7` is generated before production MD.


---

# CPPTRAJ Analysis Workflow Following My Scripts

## 22. Analysis Overview

After production MD, my analysis workflow uses `cpptraj` to:

1. Combine production trajectories
2. Strip water and ions
3. Generate a stripped topology
4. Generate a reference structure
5. Prepare a trajectory window for MM/GBSA or MM/PBSA
6. Calculate RMSD
7. Calculate RMSF
8. Calculate hydrogen bonds, if `hbond.in` is available
9. Calculate native contacts
10. Extract representative PDB frames
11. Extract the final PDB frame

The main analysis script is:

```bash
cpptraj-analysis.sh
```

It runs the analysis files in this order:

```text
sumtotal.in
summmpbsa.in
rmsd.in
rmsf.in
hbond.in
contact_5.in
run-extract.in
run-Last-PDB.in
```

> **Important:** Your master script calls `hbond.in`, but this file was not included in the uploaded analysis set. If `hbond.in` is missing in the working directory, the script will stop at Step 5 unless you create `hbond.in` or temporarily comment out that section.

---

## 23. Required Files Before Analysis

Before running the analysis pipeline, check that the following files exist:

```bash
ls addWAT.top
ls min_all.restrt
ls complex_md1.nc
ls complex_md250.nc
ls cpptraj-analysis.sh
ls sumtotal.in
ls summmpbsa.in
ls rmsd.in
ls rmsf.in
ls contact_5.in
ls run-extract.in
ls run-Last-PDB.in
```

If you only ran production up to `md50`, then files such as `complex_md226.nc` or `complex_md250.nc` will not exist. In that case, change the trajectory ranges in `sumtotal.in` and `summmpbsa.in`.

---

## 24. Master Analysis Script: `cpptraj-analysis.sh`

### Purpose

This script runs all `cpptraj` analysis steps in the correct order.

```bash
#!/bin/bash
# =============================================================================
# Master analysis script — runs all cpptraj analyses in correct order
# Usage: bash run_analysis.sh
# =============================================================================

CPPTRAJ=$AMBERHOME/bin/cpptraj

echo "============================================"
echo " MD Analysis Pipeline"
echo "============================================"

# --- Step 1: Combine all trajectories and strip water ---
echo "[1/8] Combining trajectories and stripping water (sumtotal.in)..."
$CPPTRAJ -i sumtotal.in > sumtotal.log 2>&1
echo "      Done -> sum_MD.nc, ref-complex.rst7"

# --- Step 2: Prepare MMPBSA trajectory window ---
echo "[2/8] Preparing MMPBSA trajectory window (summmpbsa.in)..."
$CPPTRAJ -i summmpbsa.in > summmpbsa.log 2>&1
echo "      Done -> mmpbsa_450_500.nc"

# --- Step 3: RMSD ---
echo "[3/8] Calculating RMSD (rmsd.in)..."
$CPPTRAJ -i rmsd.in > rmsd.log 2>&1
echo "      Done -> rmsd_complex.dat, rmsd_protein.dat, rmsd_ligand.dat"
echo "              rmsd_backbone.dat, rmsd_bindingsite.dat"

# --- Step 4: RMSF ---
echo "[4/8] Calculating RMSF (rmsf.in)..."
$CPPTRAJ -i rmsf.in > rmsf.log 2>&1
echo "      Done -> rmsf_byres.dat, rmsf_byres_Bfac.dat, rmsf_byatom_Bfac.dat"

# --- Step 5: Hydrogen bonds ---
echo "[5/8] Calculating hydrogen bonds (hbond.in)..."
$CPPTRAJ -i hbond.in > hbond.log 2>&1
echo "      Done -> nhb_all.dat, avghb_all.dat"

# --- Step 6: Native contacts ---
echo "[6/8] Calculating native contacts (contact_5.in)..."
$CPPTRAJ -i contact_5.in > contact_5.log 2>&1
echo "      Done -> contacts_5.dat, numcontacts_5.dat, contact_5.pdb"

# --- Step 7: Extract PDB frames ---
echo "[7/8] Extracting PDB frames every 10th frame (run-extract.in)..."
mkdir -p All-PDB
$CPPTRAJ -i run-extract.in > run-extract.log 2>&1
echo "      Done -> All-PDB/PDB.pdb.*"

# --- Step 8: Extract last frame ---
echo "[8/8] Extracting last frame (run-Last-PDB.in)..."
$CPPTRAJ -i run-Last-PDB.in > run-Last-PDB.log 2>&1
echo "      Done -> Last-NS.pdb"

echo "============================================"
echo " All analyses complete!"
echo " Check *.log files if any step fails."
echo "============================================"
```

### Run

```bash
bash cpptraj-analysis.sh
```

### Output Logs

Each analysis step writes a log file:

```text
sumtotal.log
summmpbsa.log
rmsd.log
rmsf.log
hbond.log
contact_5.log
run-extract.log
run-Last-PDB.log
```

If any step fails, inspect the corresponding `.log` file.

Example:

```bash
tail -n 50 rmsd.log
tail -n 50 contact_5.log
```

---

## 25. Step 1 — Combine Production Trajectories: `sumtotal.in`

### Purpose

`sumtotal.in` loads all production trajectories, calculates SASA for selected binding-site residues, strips water and ions, centers/images the system, writes the combined trajectory, and generates a stripped reference structure.

```text
parm addWAT.top

# Load all production MD trajectories (chunks 1-250)
trajin complex_md{1..250}.nc

# SASA of binding site residues
surf :124,142,143,246,247,250,253,254,264,265,266,277,278,282,290 out sasa_all.dat

# Strip water and ions, center, image, write combined trajectory
strip :WAT,Na+,Cl- outprefix stripped
center origin :1-141
image origin center familiar
autoimage familiar
trajout sum_MD.nc NetCDF
run

# Write stripped reference structure from minimization restart
clear trajin trajout
trajin min_all.restrt
strip :WAT,Na+,Cl-
trajout ref-complex.rst7
run
```

### Main Inputs

```text
addWAT.top
complex_md1.nc ... complex_md250.nc
min_all.restrt
```

### Main Outputs

| Output | Purpose |
|---|---|
| `sum_MD.nc` | Combined stripped production trajectory |
| `stripped.addWAT.top` | Topology after stripping water and ions |
| `ref-complex.rst7` | Stripped reference structure |
| `sasa_all.dat` | SASA of selected binding-site residues |

### Important Range

Your file uses:

```text
trajin complex_md{1..250}.nc
```

This assumes that production has completed from:

```text
complex_md1.nc → complex_md250.nc
```

If only `md1–md50` are available, change it to:

```text
trajin complex_md{1..50}.nc
```

### Scientific Interpretation of SASA

The `surf` command calculates solvent-accessible surface area for selected residues. In your script, the selected residues are:

```text
124,142,143,246,247,250,253,254,264,265,266,277,278,282,290
```

These appear to represent binding-site residues. SASA can help evaluate whether the binding pocket becomes more exposed or buried during simulation.

---

## 26. Step 2 — Prepare MM/GBSA or MM/PBSA Window: `summmpbsa.in`

### Purpose

`summmpbsa.in` extracts a selected production window for MM/GBSA or MM/PBSA analysis.

```text
parm addWAT.top

# Load production MD trajectories for MMPBSA window (chunks 226-250 = 450-500 ns)
trajin complex_md{226..250}.nc

# Strip water and ions, center, image, write MMPBSA trajectory
strip :WAT,Na+,Cl- outprefix stripped
center origin :1-141
image origin center familiar
autoimage familiar
trajout mmpbsa_450_500.nc NetCDF
run

# Write stripped reference restart for MMPBSA
clear trajin trajout
trajin min_all.restrt
strip :WAT,Na+,Cl-
trajout mmpbsa_450_500.rst7
run
```

### Main Inputs

```text
addWAT.top
complex_md226.nc ... complex_md250.nc
min_all.restrt
```

### Main Outputs

| Output | Purpose |
|---|---|
| `mmpbsa_450_500.nc` | Stripped trajectory window for MM/GBSA or MM/PBSA |
| `mmpbsa_450_500.rst7` | Stripped reference restart for the same system |

### Time Window Explanation

If each production chunk is 2 ns:

```text
md226–md250 = 25 chunks
25 × 2 ns = 50 ns
```

The chunk range corresponds approximately to:

```text
450–500 ns
```

### Important Warning

This file assumes that production has reached `complex_md250.nc`. If production only reached `md50`, this file will fail.

For a 100 ns simulation using `md1–md50`, you might change the range to the last 50 ns:

```text
trajin complex_md{26..50}.nc
trajout mmpbsa_50_100.nc NetCDF
trajout mmpbsa_50_100.rst7
```

---

## 27. Step 3 — RMSD Analysis: `rmsd.in`

### Purpose

`rmsd.in` calculates RMSD for:

- whole complex
- protein only
- ligand only
- backbone atoms
- selected binding-site residues

```text
parm stripped.addWAT.top
trajin sum_MD.nc

# Center, image and align trajectory
center :1-302 mass origin
image origin center
autoimage familiar

# Reference structure for RMSD calculation
reference ref-complex.rst7

# RMSD calculations — all written to separate output files for clarity
rms reference :1-302              out rmsd_complex.dat          # Whole complex
rms reference :1-301              out rmsd_protein.dat          # Protein only
rms reference :302                out rmsd_ligand.dat           # Ligand only
rms reference :1-302@CA,C,O,N     out rmsd_backbone.dat         # Backbone only
rms reference :124,142,143,246,247,250,253,254,264,265,266,277,278,282,290@CA,C,O,N  out rmsd_bindingsite.dat  # Binding site residues

go
```

### Outputs

| Output File | Meaning |
|---|---|
| `rmsd_complex.dat` | RMSD of whole complex |
| `rmsd_protein.dat` | RMSD of protein residues `1–301` |
| `rmsd_ligand.dat` | RMSD of ligand residue `302` |
| `rmsd_backbone.dat` | RMSD of backbone atoms |
| `rmsd_bindingsite.dat` | RMSD of selected binding-site residues |

### Interpretation

- **Protein RMSD** indicates overall structural stability.
- **Backbone RMSD** is usually more stable than whole-complex RMSD.
- **Ligand RMSD** indicates whether the ligand maintains a similar pose.
- **Binding-site RMSD** shows local binding-pocket stability.

### Scientific Caution

RMSD stability does not prove binding affinity. It only suggests that the selected atoms remain structurally similar to the reference over the simulated time.

---

## 28. Step 4 — RMSF Analysis: `rmsf.in`

### Purpose

`rmsf.in` calculates per-residue and per-atom fluctuations after aligning the trajectory to the reference.

```text
parm stripped.addWAT.top
trajin sum_MD.nc

reference ref-complex.rst7 [ref]

# Image and align trajectory to reference
autoimage
rms reference [ref] @C,CA,N

# RMSF per residue — backbone atoms only
atomicfluct out rmsf_byres.dat       @C,CA,N byres

# B-factor per residue (for VMD/PyMOL visualization)
atomicfluct out rmsf_byres_Bfac.dat  @C,CA,N byres   bfactor

# B-factor per atom
atomicfluct out rmsf_byatom_Bfac.dat @C,CA,N byatom  bfactor

go
```

### Outputs

| Output File | Meaning |
|---|---|
| `rmsf_byres.dat` | Per-residue backbone RMSF |
| `rmsf_byres_Bfac.dat` | Per-residue fluctuation formatted as B-factor-like output |
| `rmsf_byatom_Bfac.dat` | Per-atom fluctuation formatted as B-factor-like output |

### Interpretation

- High RMSF often indicates flexible loops, termini, or unstable regions.
- Binding-site RMSF may suggest whether ligand binding stabilizes or destabilizes local residues.
- RMSF should be interpreted with structural visualization and contact/hydrogen-bond analysis.

---

## 29. Step 5 — Hydrogen Bond Analysis: `hbond.in`

Your master script includes:

```bash
$CPPTRAJ -i hbond.in > hbond.log 2>&1
```

However, `hbond.in` was not included with the uploaded analysis files.

### Suggested Example `hbond.in`

If you want to analyze ligand-protein hydrogen bonds, create:

```text
parm stripped.addWAT.top
trajin sum_MD.nc

autoimage
reference ref-complex.rst7

# Ligand-protein hydrogen bonds
hbond LigProtHB :302 :1-301 \
  out nhb_all.dat \
  avgout avghb_all.dat \
  series

run
```

### Outputs

| Output | Meaning |
|---|---|
| `nhb_all.dat` | Number of hydrogen bonds over time |
| `avghb_all.dat` | Average hydrogen-bond information |

### Interpretation

Hydrogen-bond occupancy is more meaningful than a single-frame interaction. A stable hydrogen bond should persist across a meaningful portion of the trajectory.

---

## 30. Step 6 — Native Contact Analysis: `contact_5.in`

### Purpose

`contact_5.in` calculates native contacts between ligand residue `302` and protein residues `1–301` using a 5.0 Å cutoff and excluding hydrogens.

```text
parm stripped.addWAT.top
trajin sum_MD.nc

autoimage
reference ref-complex.rst7

# Native contacts between ligand (:302) and protein (:1-301)
# distance 5.0 A cutoff, excluding hydrogens
# writecontacts — per-contact list, out — contact count per frame

nativecontacts :302&!@H= :1-301&!@H= \
  writecontacts contacts_5.dat \
  out numcontacts_5.dat \
  reference distance 5.0 \
  mindist maxdist \
  contactpdb contact_5.pdb

run
```

### Outputs

| Output File | Meaning |
|---|---|
| `contacts_5.dat` | Per-contact list |
| `numcontacts_5.dat` | Number of contacts per frame |
| `contact_5.pdb` | Contact representation PDB |

### Interpretation

Native contacts help evaluate whether the ligand remains close to key protein residues throughout the simulation.

A high and persistent contact count may support pose stability, but it does not prove binding free energy or inhibitory activity.

---

## 31. Step 7 — Extract PDB Frames: `run-extract.in`

### Purpose

`run-extract.in` extracts every 10th frame from the combined trajectory and writes individual PDB files to the `All-PDB` folder.

```text
parm stripped.addWAT.top

# Extract every 10th frame from full trajectory as individual PDB files
# Output saved to All-PDB/ directory (must exist before running)
trajin sum_MD.nc 1 last 10

reference sum_MD.nc 1

# Center, image and align to first frame
center :1-140 mass origin
autoimage
rms reference :1-140@CA,C,N

# Write each frame as separate PDB — named PDB.pdb.1, PDB.pdb.2, ...
trajout All-PDB/PDB.pdb pdb multi

go
```

### Run

```bash
mkdir -p All-PDB
cpptraj -i run-extract.in
```

### Output

```text
All-PDB/PDB.pdb.1
All-PDB/PDB.pdb.2
All-PDB/PDB.pdb.3
...
```

### Use Case

These extracted structures can be used for:

- visual inspection
- movie generation
- representative frame analysis
- checking ligand movement
- checking structural artifacts

---

## 32. Step 8 — Extract Last Frame: `run-Last-PDB.in`

### Purpose

`run-Last-PDB.in` extracts only the last frame of `sum_MD.nc`.

```text
parm stripped.addWAT.top

# Extract only the last frame of the trajectory as a PDB file
trajin sum_MD.nc last last 1

reference sum_MD.nc 1

# Center, image and align to first frame
center :1-140 mass origin
autoimage
rms reference :1-140@CA,C,N

# Write last frame as PDB
trajout Last-NS.pdb

go
```

### Output

```text
Last-NS.pdb
```

### Use Case

`Last-NS.pdb` can be used to inspect the final structure after the full production simulation.

---

## 33. Optional Conversion: `sumMD-rms.in`

### Purpose

`sumMD-rms.in` reads `sum_MD.nc` and writes an AMBER ASCII trajectory file:

```text
sum_MD.mdcrd
```

```text
#Read in the trajectory file
# starting at 1
# 100 total snapshots
# consider each one

parm stripped.addWAT.top

trajin sum_MD.nc 1 50000 10

#reference ../../3MD_run_newest/3U1I_minALL.crd
#strip :WAT,Na+ outprefix stripped
center origin :1-302
image origin center familiar

trajout sum_MD.mdcrd 

go
```

### Output

```text
sum_MD.mdcrd
```

### When to Use

Use this only if another tool requires `.mdcrd` format. For most modern AMBER analysis, NetCDF `.nc` is preferred because it is smaller and more reliable.

---

## 34. Critical Residue-Range Consistency Check

Your analysis files use different residue ranges:

| File | Selection |
|---|---|
| `sumtotal.in` | centers `:1-141` |
| `summmpbsa.in` | centers `:1-141` |
| `rmsd.in` | complex `:1-302`, protein `:1-301`, ligand `:302` |
| `contact_5.in` | ligand `:302`, protein `:1-301` |
| `run-extract.in` | centers and aligns `:1-140` |
| `run-Last-PDB.in` | centers and aligns `:1-140` |
| `sumMD-rms.in` | centers `:1-302` |

### Why This Matters

These selections must match the actual system. Otherwise:

- RMSD may be calculated on the wrong atoms.
- Ligand residue may be misidentified.
- Protein-only analysis may include or exclude wrong residues.
- Extracted PDB frames may be aligned differently from RMSD/RMSF results.
- MM/GBSA trajectory may not match the intended complex definition.

### What to Do Before Analysis

Run:

```bash
cpptraj -p stripped.addWAT.top
```

Then inside `cpptraj`:

```text
resinfo
atominfo :302
```

Also check the solvated topology:

```bash
cpptraj -p addWAT.top
```

Then:

```text
resinfo
```

Confirm:

```text
Protein residues = ?
Ligand residue   = ?
Complex residues = ?
Binding-site residues = ?
```

After confirming, update all scripts consistently.

---

## 35. Running the Full Analysis Pipeline

### Step 1: Make the script executable

```bash
chmod +x cpptraj-analysis.sh
```

### Step 2: Run the script

```bash
bash cpptraj-analysis.sh
```

### Step 3: Check logs

```bash
tail -n 50 sumtotal.log
tail -n 50 summmpbsa.log
tail -n 50 rmsd.log
tail -n 50 rmsf.log
tail -n 50 contact_5.log
tail -n 50 run-extract.log
tail -n 50 run-Last-PDB.log
```

If `hbond.in` is missing, either create it or temporarily edit `cpptraj-analysis.sh` and comment out:

```bash
$CPPTRAJ -i hbond.in > hbond.log 2>&1
```

---

## 36. Analysis Output Summary

| Analysis | Input Script | Main Output |
|---|---|---|
| Combine trajectory | `sumtotal.in` | `sum_MD.nc`, `ref-complex.rst7`, `stripped.addWAT.top` |
| MM/GBSA window | `summmpbsa.in` | `mmpbsa_450_500.nc`, `mmpbsa_450_500.rst7` |
| RMSD | `rmsd.in` | `rmsd_complex.dat`, `rmsd_protein.dat`, `rmsd_ligand.dat`, `rmsd_backbone.dat`, `rmsd_bindingsite.dat` |
| RMSF | `rmsf.in` | `rmsf_byres.dat`, `rmsf_byres_Bfac.dat`, `rmsf_byatom_Bfac.dat` |
| Hydrogen bonds | `hbond.in` | `nhb_all.dat`, `avghb_all.dat` |
| Native contacts | `contact_5.in` | `contacts_5.dat`, `numcontacts_5.dat`, `contact_5.pdb` |
| Extract PDBs | `run-extract.in` | `All-PDB/PDB.pdb.*` |
| Extract last frame | `run-Last-PDB.in` | `Last-NS.pdb` |
| Optional mdcrd | `sumMD-rms.in` | `sum_MD.mdcrd` |

---

## 37. Method Description for CPPTRAJ Analysis

```text
Trajectory analysis was performed using cpptraj from AmberTools. Production trajectories were combined and processed by removing water molecules and ions, centering and imaging the system, and writing a combined NetCDF trajectory. A stripped reference structure was generated from the minimized restart file. RMSD was calculated for the whole complex, protein, ligand, backbone atoms, and selected binding-site residues. RMSF was calculated for backbone atoms on a per-residue basis. Ligand-protein native contacts were calculated using a 5.0 Å cutoff excluding hydrogen atoms. Representative PDB frames and the final trajectory frame were extracted for structural visualization. A selected production window was also prepared for MM/GBSA or MM/PBSA analysis.
```

---

## 38. Final Analysis Recommendation

Before trusting any analysis result, confirm these four things:

```text
1. The trajectory range in each script matches the MD chunks that actually exist.
2. The residue masks match your real system.
3. The stripped topology matches the stripped trajectory.
4. The reference structure was generated from the same topology and atom order.
```

The most important correction to check in your current analysis workflow is residue numbering. Your scripts use both `:1-141` and `:1-302` style selections. These may both be correct for different systems, but they should not be mixed accidentally.

