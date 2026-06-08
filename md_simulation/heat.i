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
 tempi=10.0,
 temp0=345.15,
 ig=-1,                          ! Randomize RNG seed
 ntr=1,                          ! Positional restraints on
 restraint_wt=5.0,               ! kcal/mol/A**2 restraint force constant
 restraintmask='!:WAT&@CA,N,C,O,ZN',
 nmropt=1,                       ! NMR restraints for temperature ramp
/
&wt
  TYPE='TEMP0', ISTEP1=0, ISTEP2=100000,
  VALUE1=10.0, VALUE2=345.15,
/
&wt TYPE='END' /
