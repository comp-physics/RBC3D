#!/bin/bash

# salloc a node before you run this because petsc configure uses srun
ml gcc/12.3.0 mvapich2/2.3.7-1 hdf5 intel-oneapi-mkl/2023.1.0 python/3.10.10 fftw/3.3.10-mva2 cmake

# # building and installing petsc 3.19.6 in packages directory
# rm -fr packages
# mkdir packages
# cd packages

# export ROOTDIR=$(pwd)
# export SRCNCDF=${ROOTDIR}/srcNETCDF
# export INSNCDF=${ROOTDIR}/NETCDF_INST

# mkdir $SRCNCDF
# mkdir $INSNCDF

# cd $SRCNCDF

# wget downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
# wget downloads.unidata.ucar.edu/netcdf-fortran/4.6.1/netcdf-fortran-4.6.1.tar.gz

# tar -xf netcdf-c-4.9.2.tar.gz
# tar -xf netcdf-fortran-4.6.1.tar.gz

# cd netcdf-c-4.9.2
# export CPPFLAGS="-DNDEBUG -DgFortran"
# export CFLAGS="-O"
# export FFLAGS="-O -w"
# export BIN=Linux2_x86_64gfort

# ./configure --prefix=${INSNCDF} --disable-netcdf-4 --disable-dap
# make all check install

# cd ../netcdf-fortran-4.6.1
# export NCDIR=${INSNCDF}
# export NFDIR=${INSNCDF}
# export CPPFLAGS=$CPPFLAGS" -I${NCDIR}/include"
# export LDFLAGS=$LDFLAGS" -L${NCDIR}/lib"
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${NCDIR}/lib
# export LIBS=$LIBS" -lnetcdf"

# # copy necessary updated files from this PR until they release new version
# # https://github.com/Unidata/netcdf-fortran/pull/429
# cd fortran
# rm -f module_netcdf_nc_data.F90
# rm -f module_typesizes.F90
# rm -f netcdf.F90
# cp ../../../../install/scripts/module_netcdf_nc_data.F90 ./
# cp ../../../../install/scripts/module_typesizes.F90 ./
# cp ../../../../install/scripts/netcdf.F90 ./
# cp ../../../../install/scripts/typesizes.mod ./
# cd ..

# ./configure --prefix=${INSNCDF}
# make check
# make install

# cd ../..

# wget https://ftp.mcs.anl.gov/pub/petsc/petsc-3.19.tar.gz
# tar -xf petsc-3.19.tar.gz

cd packages/petsc-3.19.6

# using mpiexec here instead of srun but srun works if you salloc a node (last time i (suzan) checked?)
./configure --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=mpif90 \
    --with-fortran-datatypes \
    --with-debugging=0 \
    --COPTFLAGS=-g -O3 -march=native -mtune=native \
    --CXXOPTFLAGS=-g -O3 -march=native -mtune=native \
    --FOPTFLAGS=-g -O3 -march=native -mtune=native \
    --with-blaslapack-dir=$MKLROOT \
    --with-mpiexec=mpiexec \
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