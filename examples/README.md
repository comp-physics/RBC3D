# Case Input Parameters

`case/Input/tube.in` has several input parameters that can be used to set up and control a simulation. `nCellTypes` is the number of different cells a case will support. `refRad` and `viscRat` can be set for each cell type but aren't used in the simulation. `refRad` is the equivalent spherical radius of a cell and has to be set in `initcond.F90` and `tube.F90` for RBC creation functions.

`Nt` is the total number of timesteps the simulation will run for. `Ts` is the time step size and may need to be adjusted depending on the simulation. `cell_out` controls how many timesteps are in between when cell coordinates get written to output `dat` files. 

`D/restart.LATEST.dat` is part of the restart procedure. The timesteps in between when `restart.LATEST.dat` gets generated is specified by `restart_out`. It contains the values of simulation variable at the latest timestep and can be used by `tube.F90` to continue a simulation from that timestep. The other input parameters are used by `common` functions and subroutines to control their behavior.
