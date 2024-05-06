# Velocity Fields for RBC3D

While other provided testcases simulate blood cells in a wall geometry, this example will calculate and output the velocities in a grid within an already-completed simulation.

## What Does This Do?

After building with `make`, the `field` program will generate a list of coordinates as a 3D grid inside the simulation's bounding box (without any regard to the location of the wall or the cells).
Then, it calculates velocities (using the same method to calculate blood cell motion) at each point in the list.

### How To Run?

After running your blood flow simulation, you should have a series of restart files named in the format: `restart{timestep}.dat` corresponding to the state of the simulation at the timestep, and `restart.LATEST.dat` is the most recent restart file.
Pick the restart file corresponding to the timestep for which you want the velocity field, and insert into the `./D/` directory.
Then, modify the `./Input/tube.in` configuration file to match the restart file's path; the current path is `D/restart.LATEST.dat`.


### Output File

After runnning the simulation, you will receive an output file `./D/field.csv`. This file contains the following comma-separated columns:

X | Y | Z | Vx | Vy | Vz

where, for each row, (X, Y, Z) denotes the location, and (Vx, Vy, Vz) is the velocity at that location.

#### Visualization in ParaView

You can easily visualize this velocity field in ParaView with the following steps after Importing and Parsing `field.csv`:
1. Apply the "Table to Points" Filter, selecting X=X, Y=Y, and Z=Z as your coordinates.
2. Apply the "Merge Vector Components" Filter, selecting X=Vx, Y=Vy, and Z=Vz for your vectors. Rename this to "velocity."
3. Apply the "Glyph" filter, selecting the orientation to be "velocity." Scale and color as necessary. 