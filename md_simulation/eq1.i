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
 tempi=345.15,                   ! Initial temperature in K
 temp0=345.15,                   ! Target temperature in K
 ig=-1,                          ! Randomize RNG seed
 ntr=1,                          ! Positional restraints on
 restraint_wt=5.0,               ! kcal/mol/A**2 restraint force constant
 restraintmask='!:WAT&@CA,N,C,O,ZN',
 nmropt=0,
/
