# RBC3D Build and Run Instructions

0. Use an appropriate computer
1. Ensure you have compilers and wrappers
2. Build libraries
3. Configure the `Makefile.in` in the base directory
4. Build `common/` and `examples/case/`

## Use an appropriate computer

* RBC3D __will not__ build on non-x86 hardware (like a new Mac M1+) at time of writing. Use an x86 machine (AMD or Intel processors).
* RBC3D has not been tested on WSL or Windows computers broadly. We do not recommend using WSL. Instead, use a Linux partition or a *nix-based computing cluster.
At Georgia Tech we have several, including ICE and PACE Phoenix. 

## Ensure you have compilers and wrappers

You will need `gcc`, `gfortran`, and a suitable MPI wrapper like `mvapich` (or the like).
* On PACE Phoenix you can issue `module load gcc mvapich2`.
* On COC-ICE you can issue `module load gcc/8.3.0 mvapich2/2.3.2` 

* To check for gfortran and gcc, see if `which gfortran` and `which gcc` return a path
* Similarly, to see if you can run MPI commands for later, see if `which mpicc` or `which mpif90` return a path

You will also need python3 for the PETSc install and pip modules. On Phoenix, you can module load it via `ml python/3.9.12-rkxvr6`. On Phoenix, you may also need to add the `~/.local/bin` directory to your PATH by adding this line to your `~/.bashrc`: 

```shell
export PATH="$PATH:$HOME/.local/bin"
```

## Build libraries

### MKL

* For mkl, you can `module load mkl` on the Phoenix cluster or `module load intel-oneapi-mkl/2023.1.0` on the ICE cluster.
* This will automatically set the `MKL_ROOT` environment variable necessary for `Makefile.in`
* You can check environment variables via `module show mkl`
* If your cluster does not have mkl, it's available for download and install [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?operatingsystem=linux&distributions=offline).
* If MKL is not available on your system, follow the LAPACK and BLAS step, and skip this step.
* Note that `MKL_LIB` options in Makefile.in may need to be changed depending on the version of mkl, but the [mkl link line advisor](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html#gs.9hbhed) should provide the correct link options.

### LAPACK and BLAS

* Move back into packages: `cd RBC3D/packages`
* Download the latest lapack (at time of writing `3.11`): `wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.11.tar.gz`
* Unpack it `tar -xvf v3.11.tar.gz`
* `cd lapack-3.11`
* Modify `make.inc.example` as needed 
   * `CC = mpicc` (line 9)
   * `FC = mpif90` (line 20)
* `mv make.inc.example make.inc`
* Build (this takes a few minutes): `make`
* Later, You will need the absolute path of `liblapack.a` and `librefblas.a` to configure `petsc-lite`
   * In my case this is `/storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/lapack-3.11/liblapack.a`
   * and `/storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/lapack-3.11/librefblas.a`
* Note that if you're choosing to install LAPACK/BLAS instead of MKL, you'll need to include `Makefile.lapack` in `examples/case/Makefile` instead of `Makefile.mkl` when you run the example case later.


### PETSc

* This depends on MKL or LAPACK/BLAS from above, don't attempt until one of those steps is finished
* Move back to `RBC3D/packages`
* Download PETSc: wget https://ftp.mcs.anl.gov/pub/petsc/petsc-3.19.tar.gz`
* Unpack `petsc-3.19`: `tar -xvf petsc-3.19.tar.gz`
* Copy configure script into petsc directory: `cp ../install/py_scripts/petsc_configure.py ./petsc-3.19.6`
* Descend into the directory: `cd petsc-3.19`
* Install pip configure: `pip3 install --user configure`
* If you installed MKL, run: `python3 petsc_configure.py --mkl-only`
* If you installed LAPACK/BLAS instead, run: `python3 petsc_configure.py --blas-lapack`
* Note the `--dryrun` option will show you the PETSc configure options.

* Build and Test:
```shell
make PETSC_DIR=`pwd` PETSC_ARCH=petsc_configure all
make PETSC_DIR=`pwd` PETSC_ARCH=petsc_configure check
``` 
* `make check` may return errors, but you can ignore these if the tests completed.

### FFTW

* Easiest if this is already available or can be loaded.
* PACE Phoenix has this available as a module: `module load fftw`
* `module show fftw` tells you where the library is.
* At time of writing, the FFTW library files live at `/usr/local/pace-apps/spack/packages/linux-rhel7-x86_64/gcc-10.3.0/fftw-3.3.10-dgx5szpp2x4fznqfuaoucmwieqxbgpg6/lib`
    * You will need this directory for the `Makefile.in`

### NETCDF

* Need to get to netcdff (netcdf-fortran)
* PACE Phoenix has this pre-installed
* `module load netcdf-c netcdf-cxx netcdf-fortran`
* Later you will need information about where `netcdf-fortran` is installed for `Makefile.in`. 
* Get this via `module show netcdf-fortran` and looking at the `NETCDF_FORTRANROOT`
* In my case, this is `/usr/local/pace-apps/spack/packages/linux-rhel7-x86_64/gcc-10.3.0/netcdf-fortran-4.5.4-yx5osuxluenmuvr3xnahmosfr3abeu2p/`

### Spherepack

* Navigate back to `RBC3D/packages`
* Fetch Spherepack: `git clone https://github.com/comp-physics/spherepack3.2.git`
* `cd spherepack3.2`
* Change `make.inc` if using non-GNU compilers
* Execute `make`

### makedepf90

* This is needed to run `make .depend` in both `common/` and `case/` during build.
* Sometimes available on systems by default
* PACE Phoenix doesn't seem to have it, so let's build it
* Descend into `RBC3D/packages`
* `git clone https://github.com/comp-physics/makedepf90.git`
* `cd makedepf90`
* Modify `Makefile`, including
  * Line 41: `prefix` (which should be your full path to the binary build, e.g. for me: `/storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D/packages/makedepf90`)
* Build: `make` 
* Install: `make install` 
* This will build the `makedepf90` binary in your `RBC3D/packages/makedepf90` directory

## Configure Makefile.in

You need to change the `Makefile.in` to locate all of these libraries!
This mostly just means changing the first line of `Makefile.in`:
```
WORK_DIR = /storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D
```
and the module directories
```bash
# Directories from loaded modules
FFTW_DIR = /usr/local/pace-apps/spack/packages/linux-rhel7-x86_64/gcc-10.3.0/fftw-3.3.10-dgx5szpp2x4fznqfuaoucmwieqxbgpg6
NETCDF_DIR = /usr/local/pace-apps/spack/packages/linux-rhel7-x86_64/gcc-10.3.0/netcdf-fortran-4.5.4-yx5osuxluenmuvr3xnahmosfr3abeu2p
```

* If you have an older version of `gfortran` then you will need to remove the `-fallow-argument-mismatch` flag on line 39 of `Makefile.in`

## Build

### Common

This is the main codebase.

* Descend into `RBC3D/common`
* Execute `make .depend`
* Execute `make`

### Case

This is an example case.

* Descend into `RBC3D/examples/case`
* Execute `make .depend`
* Execute `make`

## Run

In `case/` you should be able to
* `mpiexec -n 1 ./initcond`
* `mpiexec ./tube`
or substitute in your mpi runner like `srun`.

Note that on Phoenix, `srun` only works if you [salloc](https://gatech.service-now.com/home?id=kb_article_view&sysparm_article=KB0041998) a new node. Also, not specifying a node count via the `-n` flag will use all the cores on the machine which is necessary for a simulation with high cell counts or the parallel `./initcond` file in `examples/randomized_case`.

## Data and visualization

After running `srun ./tube` or the equivalent, you should see `x000*.dat`, `wall000*.dat`, and restart files.
You can load the `.dat` files into Paraview to visualize them.
