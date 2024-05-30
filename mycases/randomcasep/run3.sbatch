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
#SBATCH -t2:00:00
#SBATCH -q embers
#SBATCH --mail-user=smanasreh6@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "./D/plsnocrash9.log"
#SBATCH --open-mode=append

cd $SLURM_SUBMIT_DIR

ml gcc mvapich2 netcdf-c netcdf-cxx netcdf-fortran fftw

timeout 1h srun ./tube 
EX_CODE=$?

echo "SLURM_JOB_EXIT_CODE" 
echo $SLURM_JOB_EXIT_CODE

# move files to scratch, except for log file and restart.LATEST
cd D
cp * /storage/home/hcoda1/6/smanareh6/scratch/large_carotid/
rm r* x* w*
cp /storage/home/hcoda1/6/smanareh6/scratch/large_carotid/restart.LATEST.dat ./restart.LATEST.dat
cd ../

if [[ ${EX_CODE} == 124 ]]; then 
  scontrol requeue $SLURM_JOB_ID
fi
