#!/bin/bash
#SBATCH --job-name=tutorial_05_binette
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --output=logs/tutorial_05_binette_%j.out

cd $SLURM_SUBMIT_DIR

# Run Binette to refine and select best bins
eval "${ENV_CMD} bash snippets/05_binette.sh"