#!/bin/bash

rm -rf packages2
mkdir packages2
cd packages2

export ROOTDIR=$(pwd)
export SRCNCDF=${ROOTDIR}/srcNETCDF
export INSNCDF=${ROOTDIR}/NETCDF_INST

mkdir $SRCNCDF
mkdir $INSNCDF

cd $SRCNCDF

wget downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
wget downloads.unidata.ucar.edu/netcdf-fortran/4.6.1/netcdf-fortran-4.6.1.tar.gz

tar -xf netcdf-c-4.9.2.tar.gz
tar -xf netcdf-fortran-4.6.1.tar.gz
git clone https://github.com/Unidata/netcdf-fortran.git

cd netcdf-c-4.9.2
export CPPFLAGS="-DNDEBUG -DgFortran"
export CFLAGS="-O"
export FFLAGS="-O -w"
export BIN=Linux2_x86_64gfort

cd ../netcdf-fortran-4.6.1
export NCDIR=${INSNCDF}
export NFDIR=${INSNCDF}
export CPPFLAGS=$CPPFLAGS" -I${NCDIR}/include"
export LDFLAGS=$LDFLAGS" -L${NCDIR}/lib"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${NCDIR}/lib
export LIBS=$LIBS" -lnetcdf"

./configure --prefix=${INSNCDF}
make check
make install