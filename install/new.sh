#!/bin/bash

rm -rf packages
mkdir packages
cd packages

export ROOTDIR=$(pwd)
export SRCNCDF=${ROOTDIR}/srcNETCDF
export INSNCDF=${ROOTDIR}/NETCDF_INST

mkdir $SRCNCDF
mkdir $INSNCDF

cd $SRCNCDF

# wget downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
wget downloads.unidata.ucar.edu/netcdf-fortran/4.6.1/netcdf-fortran-4.6.1.tar.gz

# tar -xf netcdf-c-4.9.2.tar.gz
tar -xf netcdf-fortran-4.6.1.tar.gz

# cd netcdf-c-4.9.2
# export CPPFLAGS="-DNDEBUG -DgFortran"
# export CFLAGS="-O"
# export FFLAGS="-O -w"
# export BIN=Linux2_x86_64gfort

./configure --prefix=${INSNCDF} --disable-netcdf-4 --disable-dap
if (($?)); then
    echo "[install-ice.sh] Error: NETCDF-C configure failed."
    exit 1
fi
make all check install
if (($?)); then
    echo "[install-ice.sh] Error: NETCDF-C tests or install failed."
    exit 1
fi

# build and install netcdf-fortran in packages/NETCDF_INST
cd netcdf-fortran-4.6.1
# instructions for build with shared libraries
# https://docs.unidata.ucar.edu/netcdf-c/current/building_netcdf_fortran.html
export NCDIR=${NETCDF_CROOT}
export NFDIR=${INSNCDF}
export CPPFLAGS=$CPPFLAGS" -I${NCDIR}/include"
export LDFLAGS=$LDFLAGS" -L${NCDIR}/lib"

# copy necessary updated files from this PR until they release new version
# https://github.com/Unidata/netcdf-fortran/pull/429
cd fortran
rm -f module_netcdf_nc_data.F90
rm -f module_typesizes.F90
rm -f netcdf.F90
cp ../../../../install/scripts/module_netcdf_nc_data.F90 ./
cp ../../../../install/scripts/module_typesizes.F90 ./
cp ../../../../install/scripts/netcdf.F90 ./
cd ..

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