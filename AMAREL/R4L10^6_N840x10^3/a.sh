#!/bin/bash

#SBATCH --partition=main

#SBATCH --job-name=SM-ARR

#SBATCH --array=1-10%10

#SBATCH --cpus-per-task=1

#SBATCH --mem=2G

#SBATCH --time=60:00:00

#SBATCH --output=slurmTEST-%A.%a.out # stdout file

#SBATCH --error=slurmTEST-%A.%a.err  # stderr file

module purge

cd /scratch/cr718/R4L10^6_N840x10^3
 
echo "Running my BS array job:"

srun ./out   
