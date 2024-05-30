#!/bin/bash

#SBATCH --account=gts-sbryngelson3
#SBATCH -N2 --ntasks-per-node=24
#SBATCH --mem-per-cpu=2G
#SBATCH -t8:00:00
#SBATCH -q embers
#SBATCH --mail-user=smanasreh6@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "./run_logs/plsworkt.log"

cd $SLURM_SUBMIT_DIR

ml gcc mvapich2 netcdf-c netcdf-cxx netcdf-fortran fftw

# cd D
# rm -rf *
# cd ../

cd ../../common
make clean
make .depend
make

cd ../mycases/randomized_case
make clean
make .depend
make
# srun ./initcond
srun ./tube