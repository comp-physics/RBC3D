# RBC3D Build Instructions

0. Ensure you have compilers and wrappers
1. Build libraries
2. Cnfigure the `Makefile.inc` in the base directory
3. Build `common/` and `case/`

## Compilers

You will need `gcc` and `mpich` (or the like).
On PACE Phoenix you can issue `module load gcc mvapich2`.

## Build libraries

### MKL

* We need an old MKL, in this case `l_mkl_p_10.0.1.014/`
    * The tarball for this library and licenses for it are in the `packages/` directory
* Untar it via
   * `cd packages`
   * `tar -xvf l_mkl_p_10.0.1.014.tgz`

* Proceed with a user installation
    * `cd ../l_mkl_p_10.0.1.014`
    * `./install/install`
    * `2. Install the software without using RPM database (root password not required).`
    * Provide the absolute path to the license file when requested
    * There are two licenses provides, for example one is at located `install/mkl-licenses/l_GVJ78MLJ.lic`
* Move into the `mkl/interfaces/lapack95` directory
* Build the LAPACK95 here by
    * Modify the `makefile` via `FC=gfortran`
    * Execute `make`
    * Copy all module files from `mkl/interfaces/lapack95/obj/` into `RBC3D/mylib/include/`
        * This should include `mkl95_lapack.mod`, `mkl95_precision.mod`, `libmkl_lapack95.a`

### BLAS

* I use `BLAS-3.10.0`
* Download it and unpack in `RBC3D/packages/BLAS-3.10.0` (for example)
* Modify `make.inc`, like
    * `FC=mpif90`
* Then just `make`, which should create the library file `blas_LINUX.a`
* You will need the absolute path of `blas_LINUX.a` to configure `Makefile.inc`

### LAPACK

* I use `lapack-3.9.1`
* Download it and unpack in `RBC3D/packages/lapack-3.9.1` (for example)
* Modify `make.inc` as needed, `CC=mpicc` and `FC=mpif90` for example.
* `make` and (maybe?) `make install` as appropriate

### Valgrind

* Easiest if this is already available or can be loaded.
* PACE Phoenix has this available as a module: `module load valgrind`
* `module show valgrind` tells you where the library is.
* At time of writing, it is here: `/usr/local/pace-apps/manual/packages/valgrind/3.19.0/gcc-4.8.5`
    * You will need this path to build PETSc-lite

### PETSc-lite

* Depends on Valgrind, LAPACK, BLAS
* Untar (`tar -xvf <file>`) the tarbar located in `install/` to `packages/`
* Set up your environment via the environment variables like this
(You will need to put these in your shell rc, e.g., your `.bashrc`. You will not be able to build RBC3D or PETSc-lite without this)
```
export PETSC_DIR=/storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D/packages/petsc-3.0.0-p3
export PETSC_ARCH=linux-gnu-c-debug
```
* Build by first creating a new directory `mypetsc/`, for example, as `RBC3D/packages/mypetsc`
* Configure via something like this, using your own absolute paths:
```
./configure --with-cc=mpicc --with-fc=mpif90 --with-blas-lib=/storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/BLAS-3.10.0/blas_LINUX.a --with-lapack-lib=/storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/lapack-3.9.1/liblapack.a --with-valgrind-dir=/usr/local/pace-apps/manual/packages/valgrind/3.19.0/gcc-4.8.5 --prefix=/storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/mypetsc --with-shared=0 --with-mpiexec=srun --with-x11=0 --with-x=0 --with-windows-graphics=0
```
* Then `make` and `make install`

### FFTW

* Easiest if this is already available or can be loaded.
* PACE Phoenix has this available as a module: `module load fftw`
* `module show fftw` tells you where the library is.
* At time of writing, the FFTW library files live at `/usr/local/pace-apps/spack/packages/linux-rhel7-x86_64/gcc-10.3.0/fftw-3.3.10-dgx5szpp2x4fznqfuaoucmwieqxbgpg6/lib`
    * You will need this directory for the `Makefile.inc`

### LAPACK95

* Download the lapack95 interface for Fortran, for example [here](https://netlib.org/lapack95/lapack95.tgz), via `wget https://netlib.org/lapack95/lapack95.tgz`
* Untar it, using, e.g., `tar -xvf <lapack95.tgz>`
* Modify `make.inc` as needed:
    * Change this line `MODLIB = -I/storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D/mylib/include`
    * Change this one to point to your lapack-3.9.1 install:
        * `LAPACK_PATH = /storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D/packages/lapack-3.9.1`
    * Also change these `FC = gfortran -ffree-form` and `FC1 = gfortran -ffixed-form` and `OPTS0 = `
* `cd SRC`
* Execute `make single_double_complex_dcomplex`
* `cd ..`
* Execute `mv LAPACK95.a liblapack95.a`

### NETCDF

* Need to get to netcdff (netcdf-fortran)
* PACE Phoenix has this pre-installed
* `module load netcdf-c netcdf-cxx netcdf-fortran`
* You can get the information about where `netcdf-fortran` is installed, which you will need for `Makefile.inc`, via `module show netcdf-fortran`

### Spherepack

* This is located in `RBC3D/install/spherepack3.2`
* Change `make.inc` as appropriate
    * For example, `LIB=../lib/libspherepack.a`
    * And this snippet
```
ifeq ($(UNAMES),Linux)

  PGI := $(shell mpif90 2>&1)

  ifeq ($(PGI),pgf90-Warning-No files to process)

    F90 := pgf90 -module ../lib -I../lib
    CPP := pgf90 -E

  else

    F90 := gfortran -g -std=legacy
    # -fmod=../lib -I../lib
    CPP := gcc

  endif

  MAKE := gmake
  AR := /usr/bin/ar

endif
```
* Execute `make all`

### makedepf90

* This is needed to run `make .depend` in both `common/` and `case/` during build.
* Sometimes available on systems by default
* PACE Phoenix doesn't seem to have it.
* You can grab it from [here](https://github.com/outpaddling/makedepf90)
    * Modify `Makefile`, including
        * `CC=gcc` (or something to this effect)
        * `prefix` (which should be set to the directory you unpacked this thing in)
    * Then just `make` and `make install` (maybe?) should be enough to build the `makedepf90` binary in whatever directory you unpacked the source in.

## Configure Makefile.inc

You need to change the `Makefile.inc` to locate all of these libraries!
All-together, you should have something like this

```
WORK_DIR = /storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D
PETSC_DIR = /storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/mypetsc
PETSC_ARCH = linux-gnu-c-debug
include $(PETSC_DIR)/conf/variables
MKL_DIR = /storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D/packages/mkl/lib/em64t
MYLIB_DIR = /storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D/mylib

vpath $(WORK_DIR)/common

COMMON_INCLUDE = -I$(WORK_DIR)/common -I$(MYLIB_DIR)/include
INCLUDE = $(COMMON_INCLUDE) $(PETSC_INCLUDE) -I/storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D/packages/petsc-3.0.0-p3/include -I/usr/local/pace-apps/spack/packages/linux-rhel7-x86_64/gcc-10.3.0/netcdf-fortran-4.5.4-yx5osuxluenmuvr3xnahmosfr3abeu2p/include

COMMON_LIB = $(WORK_DIR)/common/libcommon.a
SPHPK_LIB = -L/storage/home/hcoda1/6/sbryngelson3/p-sbryngelson3-0/RBC3D/install/spherepack3.2/lib -lspherepack
LAPACK95_LIB = -L/storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/lapack95 -llapack95
FFTW_LIB = -L/usr/local/pace-apps/spack/packages/linux-rhel7-x86_64/gcc-10.3.0/fftw-3.3.10-dgx5szpp2x4fznqfuaoucmwieqxbgpg6/lib -lfftw3
NETCDF_LIB = -L/usr/local/pace-apps/spack/packages/linux-rhel7-x86_64/gcc-10.3.0/netcdf-fortran-4.5.4-yx5osuxluenmuvr3xnahmosfr3abeu2p/lib -lnetcdff
PETSC_LIB = -L$(PETSC_DIR)/lib $(PETSC_KSP_LIB_BASIC)
MKL_LIB = -L$(MKL_DIR) -Wl,-rpath,$(MKL_DIR) \
    -lmkl_lapack -lmkl -lguide -lpthread

# Compiler and linker
FC = mpif90
NICE = -fallow-argument-mismatch -freal-4-real-8
DEBUG =
OPTS = -O0

FFLAGS = $(NICE) $(DEBUG) $(OPTS) $(PROF) $(INCLUDE)
LDFLAGS = $(NICE) $(DEBUG) $(OPTS) $(PROF) -Wl,--copy-dt-needed-entries

%.o : %.F90
    $(FC) $(FFLAGS) -c $<
%.o : %.F
    $(FC) $(FFLAGS) -c $<

# Archiving
AR  = ar ru
RANLIB  = ranlib
```

## Build

### Common

This is the main code-base.

* Descend into `RBC3D/common`
* Modify the `Makefile` to use the path for your `makedepf90`
    * Like this: `/storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/makedepf90/makedepf90 *.F90 > .depend`
* Run `rm -f .depend` and then `make .depend`
* Execute `make`

### Case

This is an example case.

* Modify the `Makefile` again, like
```
include ../Makefile.in

LIB = $(COMMON_LIB) $(SPHPK_LIB) $(FFTW_LIB) $(PETSC_LIB) $(NETCDF_LIB) $(MKL_LIB) $(LAPACK95_LIB) $(STATIC)
EXECUTABLES = tube initcond

all : $(EXECUTABLES)

lib : $(wildcard $(WORK_DIR)/common/*.F90)
    make -C $(WORK_DIR)/common

$(EXECUTABLES) : % : %.o $(wildcard $(COMMON_DIR)/*.F90)
    make lib
    $(FC) $(LDFLAGS) -o $@ $< $(LIB)

clean :
    -rm -f $(EXECUTABLES) *.o *.mod *.a core

# Dependency
.depend : $(wildcard *.F90)
    /storage/coda1/p-sbryngelson3/0/sbryngelson3/RBC3D/packages/makedepf90/makedepf90 $^ > .depend

include .depend
```
