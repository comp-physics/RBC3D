#!/bin/bash

#run this in the case directory of RBC3D

#SBATCH --account=gts-sbryngelson3
#SBATCH -N8 --ntasks-per-node=24
#SBATCH --mem-per-cpu=4GB
#SBATCH -t8:00:00
#SBATCH -q embers
#SBATCH --mail-user=smanasreh6@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "./run_logs/carotidweb4.log"

cd $SLURM_SUBMIT_DIR

ml gcc mvapich2 netcdf-c netcdf-cxx netcdf-fortran fftw

# rm -r D/*

cd ../../common
make clean
make 

make clean
make .depend
make

cd ../examples/carotid_web
# srun ./initcond
srun ./tube
# timeout 7h srun ./tube 
# EX_CODE=$?

# echo "SLURM_JOB_EXIT_CODE" 
# echo $SLURM_JOB_EXIT_CODE

# move files to scratch, except for log file and restart.LATEST
# cp -r D ~/../../../../scratch1/4/smanasreh6/carotid_web2
# cd D
# rm restart0* x*
# cd ../

# sbatch run.sbatch

# if [[ ${EX_CODE} == 124 ]]; then 
#   sbatch run.sbatch
#   # scontrol requeue $SLURM_JOB_ID
# fi