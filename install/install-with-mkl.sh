#!/bin/bash

# salloc a node before you run this because petsc configure uses srun

# building and installing petsc 3.19.6 in packages directory
mkdir packages
cd packages

wget https://ftp.mcs.anl.gov/pub/petsc/petsc-3.19.tar.gz
tar -xf petsc-3.19.tar.gz

cp ../install/scripts/petsc_configure.py ./petsc-3.19.6
cd petsc-3.19.6

# if these configure options don't work, it's probably a path issue
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

make PETSC_DIR=`pwd` PETSC_ARCH=arch-linux-c-opt all
make PETSC_DIR=`pwd` PETSC_ARCH=arch-linux-c-opt check

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

echo "Done installing RBC3D!"