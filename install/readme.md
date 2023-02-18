# RBC3D Build Instructions

0. Ensure you have compilers and wrappers
1. Build libraries
2. Cnfigure the `Makefile.inc` in the base directory
3. Build `common/` and `case/`

## Compilers

You will need `gcc` and `mpich` (or the like).
* On PACE Phoenix you can issue `module load gcc mvapich2`.
* On COC-ICE you can issue `module load gcc/8.3.0 mvapich2/2.3.2`

## Build libraries

### MKL

* We need an old MKL, in this case `l_mkl_p_10.0.1.014/`
* Fetch it: 
   * `cd packages/`
   * `wget https://www.dropbox.com/s/gljrvl6p2f5x3go/l_mkl_p_10.0.1.014.tgz`
   * `tar -xvf l_mkl_p_10.0.1.014.tgz`
* Proceed with a user installation
    * `mkdir RBC3D/packages/mkl` 
    * `cd l_mkl_p_10.0.1.014`
    * `./install.sh`
    * `3. Install as current user.`
    * `1. Install`
    * `2. Provide the absolute path for an existing license file.`
    * `Please type a selection or License file name or port@hostname: /storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/mkl-licenses/l_GVJ78MLJ.lic`
    * `2. Install the software without using RPM database (root password not required).`
    * `Enter`
    * `Enter`
    * `d`, `d`, `d`, `d`, type `accept` and return
    * `Where do you want to install to?`: 
      * Via another terminal move into the `RBC3D/packages/mkl/` directory and issue `pwd`
      * Use this location for the prompt 
      * Accept warnings about overwriting the existing directory
    * `Enter to continue`
    * `Enter to continue`
    * Done
* Build the MKL LAPACK95 by
    * Move into the MKL LAPACK95 directory: `cd RBC3D/packages/mkl/interfaces/lapack95` 
    * Modify the `makefile` to `FC = gfortran` on lines 47 and 50.
    * Execute `make libem64t`
    * Execute `make lib64`

### BLAS

* Move back into packages: `cd RBC3D/packages`
* Download the latest BLAS (at time of writing `3.11.0`): `wget http://www.netlib.org/blas/blas-3.11.0.tgz`
* Unpack it: `tar -xvf blas-3.11.0.tgz`
* `cd BLAS-3.11.0`
* Modify `make.inc` line 18 as `FC = mpif90`
* Execute `make`, which will create the library file `blas_LINUX.a`
* Later, You will need the absolute path of `blas_LINUX.a` to configure `petsc-lite`
   * In my case this is `/storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/BLAS-3.11.0/blas_LINUX.a`

### LAPACK

* Move back into packages: `cd RBC3D/packages`
* Download the latest lapack (at time of writing `3.11`): `wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.11.tar.gz`
* Unpack it `tar -xvf v3.11.tar.gz`
* `cd lapack-3.11`
* Modify `make.inc.example` as needed 
   * `CC = mpicc` (line 9)
   * `FC = mpif90` (line 20)
* `mv make.inc.example make.inc`
* Build (this takes a few minutes): `make`
* Later, You will need the absolute path of `liblapack.a` to configure `petsc-lite`
   * In my case this is `/storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/lapack-3.11/liblapack.a`

### Valgrind

* Easiest if this is already available or can be loaded.
* PACE Phoenix has this available as a module: `module load valgrind`
* `module show valgrind` tells you where the library is.
* At time of writing, it is here: `/usr/local/pace-apps/manual/packages/valgrind/3.19.0/gcc-4.8.5`
    * You will need this path to build PETSc-lite

### PETSc-lite

* This depends on Valgrind, LAPACK, BLAS from above, don't attempt until those steps are finished
* Move back to `RBC3D/packages`
* Unpack `petsc-lite`: `tar -xvf petsc-lite-3.0.0-p3.tar.gz`
* Set up your environment via the environment variables 
   * Get the absolute path of the unpacked petsc via
      * `cd petsc-3.0.0-p3`
      * `pwd`: In my case: `/storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D/packages/petsc-3.0.0-p3`
   * If you are using something other than bash, look up how to set environment variables for it, otherwise:
      * Execute (notice this is the path from above): `export PETSC_DIR=/storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D/packages/petsc-3.0.0-p3`
      * Execute: `export PETSC_ARCH=linux-gnu-c-debug`
* Ascend up a directory and create a new build directory like `mkdir RBC3D/packages/mypetsc`
* Configure via something like this, using your own absolute paths (for blas, lapack, valgrind, and mypetsc):
```
./configure --with-cc=mpicc --with-fc=mpif90 --with-blas-lib=/storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/BLAS-3.11.0/blas_LINUX.a --with-lapack-lib=/storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/lapack-3.11/liblapack.a --with-valgrind-dir=/usr/local/pace-apps/manual/packages/valgrind/3.19.0/gcc-4.8.5 --prefix=/storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/mypetsc --with-shared=0 --with-mpiexec=srun --with-x11=0 --with-x=0 --with-windows-graphics=0
```
* Build: `make` 
* Install: `make install`

### FFTW

* Easiest if this is already available or can be loaded.
* PACE Phoenix has this available as a module: `module load fftw`
* `module show fftw` tells you where the library is.
* At time of writing, the FFTW library files live at `/usr/local/pace-apps/spack/packages/linux-rhel7-x86_64/gcc-10.3.0/fftw-3.3.10-dgx5szpp2x4fznqfuaoucmwieqxbgpg6/lib`
    * You will need this directory for the `Makefile.inc`

### LAPACK95

* Navigate back to `RBC3D/packages`
* Fetch lapack95: `wget https://netlib.org/lapack95/lapack95.tgz`
* Expand it: `tar -xvf lapack95.tgz`
* cd `LAPACK95`
* Modify `make.inc`:
    * Change line 23 of `make.inc` to point to the full path of your `lapack-3.11` build:
        * `LAPACK_PATH = /storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D/packages/lapack-3.11`
    * Also change 
      * Line 6 to: `FC = gfortran -ffree-form` 
      * Line 7 to: `FC1 = gfortran -ffixed-form` 
      * Line 16 to `OPTS0 = `
* `cd SRC`
* Execute `make single_double_complex_dcomplex`
* `cd ..` (back to `RBC3D/packages/LAPACK95`)
* Execute `mv lapack95.a liblapack95.a`

### NETCDF

* Need to get to netcdff (netcdf-fortran)
* PACE Phoenix has this pre-installed
* `module load netcdf-c netcdf-cxx netcdf-fortran`
* Later you will need information about where `netcdf-fortran` is installed for `Makefile.inc`. 
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
* `git clone https://github.com/outpaddling/makedepf90.git`
* `cd makedepf90`
* Modify `Makefile`, including
  * Line 30: `CC=gcc` (or something to this effect)
  * Line 41: `prefix` (which should be your full path to the binary build, e.g. for me: `/storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D/packages/makedepf90`)
* Build: `make` 
* Install: `make install` 
* This will build the `makedepf90` binary in your `RBC3D/packages/makedepf90` directory

## Configure Makefile.inc

You need to change the `Makefile.inc` to locate all of these libraries!
This mostly just means changing the first 14 lines of `Makefile.inc`:

```
WORK_DIR = /storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D
PETSC_DIR = $(WORK_DIR)/packages/mypetsc
include $(PETSC_DIR)/conf/variables
MKL_DIR = $(WORK_DIR)/packages/mkl
LAPACK95_DIR = $(WORK_DIR)/packages/LAPACK95
SPHEREPACK_DIR = $(WORK_DIR)/packages/spherepack3.2

# Makedependf90 binary
MAKEDEPEND_BIN = $(WORK_DIR)/packages/makedepf90/makedepf90

# Directories from loaded modules
FFTW_DIR = /usr/local/pace-apps/spack/packages/linux-rhel7-x86_64/gcc-10.3.0/fftw-3.3.10-dgx5szpp2x4fznqfuaoucmwieqxbgpg6
NETCDF_DIR = /usr/local/pace-apps/spack/packages/linux-rhel7-x86_64/gcc-10.3.0/netcdf-fortran-4.5.4-yx5osuxluenmuvr3xnahmosfr3abeu2p
```

* If you have an older version of `gfortran` then you will need to remove the `-fallow-argument-mismatch` flag on line 38 of `Makefile.in`

## Build

### Common

This is the main codebase.

* Descend into `RBC3D/common`
* Execute `make .depend`
* Execute `make`

### Case

This is an example case.

* Descend into `RBC3D/case`
* Execute `make .depend`
* Execute `make`

