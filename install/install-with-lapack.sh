#!/bin/bash

# salloc a node before you run this because petsc tests use srun
set -e -x

# # create packages directory
mkdir packages
cd packages

# build and install lapack and blas
wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.11.tar.gz
tar -xf v3.11.tar.gz
cd lapack-3.11
cp ../../install/scripts/make.inc ./
echo "PWD: `PWD`"
make -j 8

cd ..
# echo "PWD: `PWD`"
# build and install petsc 3.19.6 in packages directory
wget https://ftp.mcs.anl.gov/pub/petsc/petsc-3.19.tar.gz
tar -xf petsc-3.19.tar.gz

parentdir="$(dirname `pwd`)"
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

make PETSC_DIR=`pwd` PETSC_ARCH=petsc_configure all
make PETSC_DIR=`pwd` PETSC_ARCH=petsc_configure check

# build and install spherepack
cd ..
git clone https://github.com/comp-physics/spherepack3.2.git
cd spherepack3.2
make -j 8

# build and install makedepf90
cd ..
git clone https://github.com/comp-physics/makedepf90.git
cd ../install/scripts
python3 mdf90_replace.py
cd ../../packages/makedepf90
make -j 8
make install

echo "Done installing!"