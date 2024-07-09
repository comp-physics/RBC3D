#!/bin/bash

module load gcc/12.1.0-qgxpzk mvapich2/2.3.7-733lcv mkl python/3.9.12-rkxvr6 netcdf-fortran cmake

# create packages directory
mkdir packages
cd packages

# build and install fftw
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar -xf fftw-3.3.10.tar.gz
cd fftw-3.3.10
export FFTW_DIR=`pwd`/build/
./configure --prefix="$FFTW_DIR"
make clean
make -j
make install
cd ..

# build and install lapack and blas
wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.11.tar.gz
tar -xf v3.11.tar.gz
cd lapack-3.11
cp ../../install/scripts/make.inc ./
echo "PWD: `PWD`"
make -j 8

cd ..
# build and install petsc 3.19.6 in packages directory
wget https://ftp.mcs.anl.gov/pub/petsc/petsc-3.19.tar.gz
tar -xf petsc-3.19.tar.gz

parentdir=$(pwd)
echo "parentdir: $parentdir"

cd petsc-3.19.6
./configure --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=mpif90 \
    --with-fortran-datatypes \
    --with-debugging=0 \
    --COPTFLAGS=-g -O3 -march=native -mtune=native \
    --CXXOPTFLAGS=-g -O3 -march=native -mtune=native \
    --FOPTFLAGS=-g -O3 -march=native -mtune=native \
    --with-blas-lib=$parentdir/packages/lapack-3.11/librefblas.a \
    --with-lapack-lib=$parentdir/packages/lapack-3.11/liblapack.a \
    --with-mpiexec=srun \
    --with-shared-libraries=0 \
    --with-x11=0 --with-x=0 --with-windows-graphics=0

if (($?)); then
    echo "[install.sh] Error: PETSc configure failed. See configure.log for more details"
    cat configure.log
    exit 1
fi

make PETSC_DIR=`pwd` PETSC_ARCH=arch-linux-c-opt all
make PETSC_DIR=`pwd` PETSC_ARCH=arch-linux-c-opt check

# build and install spherepack
cd ..
git clone https://github.com/comp-physics/spherepack3.2.git
cd spherepack3.2
make -j 8

# build and install makedepf90
cd ../..
./rbc.sh install-makedepf90

echo "Done installing RBC3D!"
