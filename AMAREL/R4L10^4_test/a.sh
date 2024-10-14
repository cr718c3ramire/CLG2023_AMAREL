#!/bin/bash

#SBATCH --partition=main

#SBATCH --job-name=SM-ARR

#SBATCH --array=1-10%10

#SBATCH --cpus-per-task=1

#SBATCH --mem=2G

#SBATCH --time=48:00:00

#SBATCH --output=R4SlurmN8120-%A.%a.out # stdout file

#SBATCH --error=R4SlurmN8120-%A.%a.err  # stderr file

module purge

cd /scratch/cr718/R4L10^4_test
 
echo "Running my BS array job:"

srun ./out   
