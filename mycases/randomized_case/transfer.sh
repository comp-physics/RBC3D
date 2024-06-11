#!/bin/bash

#SBATCH --account=gts-sbryngelson3
#SBATCH -N8 --ntasks-per-node=24
#SBATCH --mem-per-cpu=2G
#SBATCH -t8:00:00
#SBATCH -q embers
#SBATCH --mail-user=smanasreh6@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "./run_logs/nocrashpls3.log"

cd $SLURM_SUBMIT_DIR

for ((i=1; i<=20; i=i+1))
do 
    sleep 30m
    rsync -a D/x0*.dat ~/scratch/bigwbcs
    rsync -a D/restart0*.dat ~/scratch/bigwbcs
done