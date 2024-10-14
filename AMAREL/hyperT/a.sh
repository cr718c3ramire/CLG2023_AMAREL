#!/bin/bash

#SBATCH --partition=main

#SBATCH --job-name=SM-ARR

#SBATCH --array=1-10%10

#SBATCH --cpus-per-task=1

#SBATCH --mem=2G

#SBATCH --time=24:00:00

#SBATCH --output=hR4L10^4N300-%A.%a.out # stdout file

#SBATCH --error=hR4L10^4N300-%A.%a.err  # stderr file

module purge

cd /scratch/cr718/hyperT
 
echo "Running my BS array job:"

srun ./out  # let's see what this does, maybe the below commented is better

#srun /scratch/cr718/./out #> $SLURM_ARRAY_TASK_ID.output
 
