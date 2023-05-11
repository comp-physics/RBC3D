#!/bin/bash

#Do the following before running
#change initcond.F90 layers variable
#change the #cpus and #cores here
#change log file name to reflect proc and layer count
#change num processes in srun to reflect the total execution threads

#run this in the case directory of RBC3D

#SBATCH --account=gts-coc-coda20
#SBATCH -N1 --ntasks-per-node=24
#SBATCH --mem-per-cpu=1G
#SBATCH -t2:00:00
#SBATCH -qinferno
#SBATCH --mail-user=smullangi3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "./run_logs/250It_RBCSim_16L_1CPU.log"

cd $SLURM_SUBMIT_DIR

cd D
rm -rf x*
rm -rf r*
rm -rf w*
cd ../

ml gcc mvapich2 netcdf-c netcdf-cxx netcdf-fortran fftw

make clean
make
srun -n 1 ./initcond
time srun -n 24 ./tube