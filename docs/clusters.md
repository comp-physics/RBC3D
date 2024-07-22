# Install on a GT Cluster

Georgia Tech has several computing clusters and using a cluster is necessary to run simulations with many cells. We've provided install scripts that can install all necessary packages on these clusters. They're used inside the github runners (linked below), so they're guaranteed to at least compile the codebase if similar instructions are followed. The same instructions can be used for other clusters, but modules will be named differently.


<p align="left">
  <a href="https://github.com/comp-physics/RBC3D/actions/workflows/phoenix.yml">
    <img src="https://github.com/comp-physics/RBC3D/actions/workflows/phoenix.yml/badge.svg" />
  </a>
  <a href="https://github.com/comp-physics/RBC3D/actions/workflows/ice.yml">
    <img src="https://github.com/comp-physics/RBC3D/actions/workflows/ice.yml/badge.svg" />
  </a>
</p>

## PACE Phoenix

To install on PACE Phoenix, you need to salloc a node to make sure `srun` is available and then run this in the RBC3D root directory: 

```shell
module load gcc/12.1.0-qgxpzk mvapich2/2.3.7-733lcv python/3.9.12-rkxvr6 netcdf-fortran cmake
./rbc.sh install-phoenix
```

Note that if the `gcc`, `mvapich2`, `mkl`, and `fftw` modules work on your Phoenix account, you should use this installer script for a faster build. You should try this one first before the other one, but it is not guaranteed to work.

```shell
module load gcc mvapich2 mkl python/3.9.12-rkxvr6 netcdf-fortran fftw cmake
./rbc.sh install
```
## PACE ICE

If you're on the ICE cluster, you can use this installer script.

```shell
module load gcc/12.3.0 mvapich2/2.3.7-1 netcdf-c hdf5/1.14.1-2-mva2 intel-oneapi-mkl/2023.1.0 python/3.10.10 fftw/3.3.10-mva2 cmake
./rbc.sh install-ice
```

## Environment Variables

Before you can run cmake, you must set `PETSC_DIR` and `PETSC_ARCH` environment variables. You can place them in your `~/.bashrc`. This path depends on where you placed RBC3D. To get the path to where you placed it you can run this from your RBC3D root directory:

```shell
rootdir=`pwd`
echo -e "export PETSC_DIR=$rootdir/packages/petsc-3.21.3 \nexport PETSC_ARCH=arch-linux-c-opt" >> ~/.bashrc
```

## Running an Example Case

Then to execute and run a case, you can:
```shell
mkdir build
cd build
cmake ..
make -j 8 case # or just `make` to make common and all the cases
cd case
srun -n 1 ./initcond
srun ./tube # command to use all the nodes
```

This will generate output files in `build/case/D`. To keep output files in `examples/case/D` and use input files in `examples/case/Input`, you can do this instead once files are built. I recommend this way.

```shell
cd examples/case
srun -n 1 ../../build/case/initcond
srun ../../build/case/tube
```
