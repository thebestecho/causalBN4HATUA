#!/bin/bash

#SBATCH --job-name="no o"
#SBATCH --partition=medium
#SBATCH --output=arrayJob_%A.out
#SBATCH --array=1001-1003
#SBATCH --cpus-per-task=1
#SBATCH --mem=3GB
#SBATCH --mail-user=xk5@st-andrews.ac.uk
#SBATCH --mail-type=END

echo "Starting job on $HOSTNAME"
echo [`date +"%Y-%m-%d %H:%M:%S"`]

Rscript --vanilla overall.R $SLURM_ARRAY_TASK_ID

echo "Job finished"
