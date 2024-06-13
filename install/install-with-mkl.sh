#!/bin/bash

# salloc a node before you run this because petsc configure uses srun

ml python/3.9.12-rkxvr6 gcc mvapich2 mkl netcdf-c netcdf-cxx netcdf-fortran fftw

# building and installing petsc 3.19.6 in packages directory
mkdir packages
cd packages

wget https://ftp.mcs.anl.gov/pub/petsc/petsc-3.19.tar.gz
tar -xf petsc-3.19.tar.gz

# echo "BEFORE pip3 install --user configure"
# pip3 install --user configure
# echo "AFTER pip3 install --user configure"

cp ../install/scripts/petsc_configure.py ./petsc-3.19.6
cd petsc-3.19.6

./configure --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=mpif90 \
    --with-fortran-datatypes \
    --with-debugging=0 \
    --COPTFLAGS=-g -O3 -march=native -mtune=native \
    --CXXOPTFLAGS=-g -O3 -march=native -mtune=native \
    --FOPTFLAGS=-g -O3 -march=native -mtune=native \
    --with-blaslapack-dir=$MKLROOT \
    --with-mpiexec=srun \
    --with-x11=0 --with-x=0 --with-windows-graphics=0
# python3 petsc_configure.py --mkl-only --dryrun
# python3 petsc_configure.py --mkl-only

make PETSC_DIR=`pwd` PETSC_ARCH=arch-linux-c-opt all
make PETSC_DIR=`pwd` PETSC_ARCH=arch-linux-c-opt check

# build and install spherepack
cd ..
git clone https://github.com/comp-physics/spherepack3.2.git
cd spherepack3.2
make

# build and install makedepf90
cd ..
git clone https://github.com/comp-physics/makedepf90.git
cd ../install/scripts
echo "PWD $(pwd)"
python3 mdf90_replace.py
cd ../../packages/makedepf90
make
make install

echo "Done installing RBC3D!"