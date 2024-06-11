#!/bin/bash

#SBATCH --account=gts-sbryngelson3
#SBATCH -N4 --ntasks-per-node=24
#SBATCH --mem-per-cpu=2G
#SBATCH -t1:00:00
#SBATCH -q embers
#SBATCH --mail-user=smanasreh6@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "./run_logs/singint.log"

cd $SLURM_SUBMIT_DIR

cd D
rm -rf *x*
rm -rf r*
rm -rf w*
cd ../

ml gcc mvapich2 mkl netcdf-c netcdf-cxx netcdf-fortran fftw

cd ../../common
make clean
make .depend
make

cd ../testcases/case
make clean
make .depend
make
srun -n 1 ./initcond
srun ./tube
