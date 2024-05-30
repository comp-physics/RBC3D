#!/bin/bash

#Do the following before running
#change initcond.F90 layers variable
#change the #cpus and #cores here
#change log file name to reflect proc and layer count
#change num processes in srun to reflect the total execution threads

#run this in the case directory of RBC3D

#SBATCH --account=gts-sbryngelson3
#SBATCH -N4 --ntasks-per-node=24
#SBATCH --mem-per-cpu=4GB
#SBATCH -t8:00:00
#SBATCH -q embers
#SBATCH --mail-user=smanasreh6@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "./run_logs/bigcasewbcs.log"

cd $SLURM_SUBMIT_DIR

cd D
rm *
cd ../

echo "NEW SBATCH FILE"

ml gcc mvapich2 netcdf-c netcdf-cxx netcdf-fortran fftw

cd ../../common
make clean
make .depend
make

cd ../mycases/randomized_case
make clean
make .depend
make
srun ./initcond #for continue run, comment this out

sbatch run.sbatch