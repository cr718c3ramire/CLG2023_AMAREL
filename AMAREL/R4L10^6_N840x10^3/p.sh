
#!/bin/bash

#SBATCH --partition=main          # Partition (job queue)

#SBATCH --requeue                 # Return job to the queue if preempted

#SBATCH --job-name=cr718test	  # Assign a short name to your job

#SBATCH --nodes=1                 # Number of nodes you require

#SBATCH --ntasks=1                # Total # of tasks across all nodes

#SBATCH --cpus-per-task=1         # Cores per task (>1 if multithread tasks)

#SBATCH --mem=16000                # Real memory (RAM) required (MB)

#SBATCH --time=02:00:00           # Total run time limit (HH:MM:SS)


cd /scratch/cr718/R4L10^5

module purge


srun ./out


