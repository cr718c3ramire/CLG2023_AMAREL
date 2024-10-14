#!/bin/bash

#SBATCH --partition=main

#SBATCH --job-name=SM-ARR

#SBATCH --array=1-10%10

#SBATCH --cpus-per-task=1

#SBATCH --mem=2G

#SBATCH --time=24:00:00

#SBATCH --output=slurmTEST-%A.%a.out # stdout file

#SBATCH --error=slurmTEST-%A.%a.err  # stderr file

module purge

cd /scratch/cr718/CLG2023

mv out amarel_runs #<---assumes out is name of compiled/executable
cd amarel_runs #<-- moved a.sh here and will run sbatch a.sh from here

echo "Running my BS array job:"

srun ./out  # let's see what this does, maybe the below commented is better

#srun /scratch/cr718/./out #> $SLURM_ARRAY_TASK_ID.output
 
