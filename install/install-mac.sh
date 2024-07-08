#!/bin/bash

mkdir packages
cd packages

# build and install fftw
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar -xf fftw-3.3.10.tar.gz
cd fftw-3.3.10
export FFTW_DIR=`pwd`/build/
./configure --prefix="$FFTW_DIR"
if (($?)); then
    echo "[install-mac.sh] Error: FFTW configure failed."
    exit 1
fi
make clean
make -j
if (($?)); then
    echo "[install-mac.sh] Error: FFTW make failed."
    exit 1
fi
make install
if (($?)); then
    echo "[install-mac.sh] Error: FFTW install failed."
    exit 1
fi
cd ..

# build and install lapack and blas
wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.11.tar.gz
tar -xf v3.11.tar.gz
cd lapack-3.11
cp ../../install/scripts/make.inc ./
make -j 24
if (($?)); then
    echo "[install-mac.sh] Error: LAPACK make failed."
    exit 1
fi
cd ..

# build and install netcdf-c in packages/NETCDF_INST
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

cd netcdf-c-4.9.2

./configure --prefix=${INSNCDF} --disable-netcdf-4 --disable-dap
if (($?)); then
    echo "[install-mac.sh] Error: NETCDF-C configure failed."
    exit 1
fi
make all check install
if (($?)); then
    echo "[install-mac.sh] Error: NETCDF-C tests or install failed."
    exit 1
fi

# build and install netcdf-fortran in packages/NETCDF_INST
cd ../netcdf-fortran-4.6.1
# instructions for build with shared libraries
# https://docs.unidata.ucar.edu/netcdf-c/current/building_netcdf_fortran.html
export NCDIR=${INSNCDF}
export NFDIR=${INSNCDF}
export CPPFLAGS=$CPPFLAGS" -I${NCDIR}/include"
export LDFLAGS=$LDFLAGS" -L${NCDIR}/lib"

cd ../../../install/scripts
python3 netcdf_replace.py

cd ../../packages/srcNETCDF/netcdf-fortran-4.6.1
cd fortran
cp typeSizes.F90 module_typesizes.F90
cd ..
# rm -f module_netcdf_nc_data.F90
# rm -f module_typesizes.F90
# rm -f netcdf.F90
# cp ../../../../install/scripts/module_netcdf_nc_data.F90 ./
# cp ../../../../install/scripts/module_typesizes.F90 ./
# cp ../../../../install/scripts/netcdf.F90 ./
# cd ..

./configure --prefix=${INSNCDF}
if (($?)); then
    echo "[install-mac.sh] Error: NETCDF-Fortran configure failed."
    exit 1
fi
make check
if (($?)); then
    echo "[install-mac.sh] Error: NETCDF-Fortran tests failed."
    exit 1
fi
make install

cd ../..

# building and installing petsc 3.21.3 in packages directory
wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.21.3.tar.gz

tar -xf petsc-3.21.3.tar.gz

packagesdir=$(pwd)
echo "packagesdir: $packagesdir"

cd petsc-3.21.3
./configure --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=mpif90 \
    --with-fortran-datatypes \
    --with-debugging=0 \
    --COPTFLAGS=-g -O3 -march=native -mtune=native \
    --CXXOPTFLAGS=-g -O3 -march=native -mtune=native \
    --FOPTFLAGS=-g -O3 -march=native -mtune=native \
    --with-blas-lib=$packagesdir/lapack-3.11/librefblas.a \
    --with-lapack-lib=$packagesdir/lapack-3.11/liblapack.a \
    --with-mpiexec=mpiexec \
    --with-x11=0 --with-x=0 --with-windows-graphics=0

if (($?)); then
    echo "[install-mac.sh] Error: PETSc configure failed. See configure.log for more details"
    cat configure.log
    exit 1
fi

make PETSC_DIR=`pwd` PETSC_ARCH=arch-darwin-c-opt all
make PETSC_DIR=`pwd` PETSC_ARCH=arch-darwin-c-opt check

# build and install spherepack
cd ..
git clone https://github.com/comp-physics/spherepack3.2.git
cd spherepack3.2
make -j 8

echo "Done installing RBC3D!"
