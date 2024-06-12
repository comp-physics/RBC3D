#!/bin/bash

#SBATCH --account=gts-sbryngelson3
#SBATCH -N4 --ntasks-per-node=24
#SBATCH --mem-per-cpu=2G
#SBATCH -t1:00:00
#SBATCH -q embers
#SBATCH --mail-user=smanasreh6@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL

cd $SLURM_SUBMIT_DIR
echo "Running in $(pwd):"

ml gcc mvapich2 mkl petsc netcdf-c netcdf-cxx netcdf-fortran fftw

mkdir packages
cd packages

# build and install spherepack
cd ..
git clone https://github.com/comp-physics/spherepack3.2.git
cd spherepack3.2
make

# build and install makedepf90
cd ..
git clone https://github.com/comp-physics/makedepf90.git
cd makedepf90
# it works without setting prefix
make
make install

cd ../../common
make .depend
make

cd ..examples/case
make .depend
make

srun -n 1 ./initcond
srun ./tube