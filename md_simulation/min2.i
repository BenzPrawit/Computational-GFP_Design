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
