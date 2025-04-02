#!/bin/bash
#SBATCH --job-name="overall networks"
#SBATCH --partition=short
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err
#SBATCH --array=1-1000
#SBATCH --mem=1GB
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=xk5@st-andrews.ac.uk
#SBATCH --mail-type=END

echo "Starting job on $HOSTNAME"
echo [`date +"%Y-%m-%d %H:%M:%S"`]

# Define an array of variables to intervene on
INTERVENED_VARIABLES=("age" "education" "workingstatus" "deprived_assets" "overcrowded" "toilet" "protect_drinkingwater" "aware_AMR")

# Loop over each variable
for INTERVENED_VARIABLE in "${INTERVENED_VARIABLES[@]}"
do
  # Execute Python script for each file index and each variable
  srun python process_network.py "$SLURM_ARRAY_TASK_ID" "$INTERVENED_VARIABLE"
done

echo "Job finished"


