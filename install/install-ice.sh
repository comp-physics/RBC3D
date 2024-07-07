#!/bin/bash

ml gcc/12.3.0 mvapich2/2.3.7-1 netcdf-c hdf5/1.14.1-2-mva2 intel-oneapi-mkl/2023.1.0 python/3.10.10 fftw/3.3.10-mva2 cmake

# build and install netcdf-c in packages/NETCDF_INST
rm -fr packages
mkdir packages
cd packages

export ROOTDIR=$(pwd)
export SRCNCDF=${ROOTDIR}/srcNETCDF
export INSNCDF=${ROOTDIR}/NETCDF_INST

mkdir $SRCNCDF
mkdir $INSNCDF

cd $SRCNCDF

wget downloads.unidata.ucar.edu/netcdf-fortran/4.6.1/netcdf-fortran-4.6.1.tar.gz

tar -xf netcdf-fortran-4.6.1.tar.gz

# build and install netcdf-fortran in packages/NETCDF_INST
cd netcdf-fortran-4.6.1
# instructions for build with shared libraries
# https://docs.unidata.ucar.edu/netcdf-c/current/building_netcdf_fortran.html
export NCDIR=${NETCDF_CROOT}
export NFDIR=${INSNCDF}
export CPPFLAGS=$CPPFLAGS" -I${NCDIR}/include"
export LDFLAGS=$LDFLAGS" -L${NCDIR}/lib"

./configure --prefix=${INSNCDF}
if (($?)); then
    echo "[install-ice.sh] Error: NETCDF-Fortran configure failed."
    exit 1
fi
make check
if (($?)); then
    echo "[install-ice.sh] Error: NETCDF-Fortran tests failed."
    exit 1
fi
make install

cd ../..

# building and installing petsc 3.19.6 in packages directory
wget https://ftp.mcs.anl.gov/pub/petsc/petsc-3.19.tar.gz
tar -xf petsc-3.19.tar.gz

cd petsc-3.19.6

# using mpiexec here instead of srun but srun works too --with-mpiexec=mpiexec \
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

echo "Done installing RBC3D!"
