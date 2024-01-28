# RBC3D Simulation Sample Files

This directory contains a collection of various samples of files which might help with creating various types of simulations.

## Sample Walls
The `./sample_walls/` subdirectory contains mesh data for various simple blood vessel geometries, which can be used to generate initial conditions for RBC3D.
Examples of how these walls are used can be seen in various cases within the `/examples/` directory. 
The `./sample_walls/README.md` also contains extra information on each wall-mesh, as well as directions on how to create a custom wall-mesh for RBC3D simulations.

## Sample Cells
The `./sample_cells/` subdirectory contains mesh data for different types of blood cells for simulation. 
These blood cells are stored as meshes and imported into the simulation, rather than only mathematically generated. 
The Fortran code within `/examples/case_diff_celltypes/` is a good example of how to use some of the cell-mesh files in a simulation. 
For more information regarding the specific sample cells, read `./sample_cells/README.md`.