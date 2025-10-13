#!/bin/bash
#SBATCH --job-name=tutorial_02_assembly
#SBATCH --cpus-per-task=12
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --output=logs/tutorial_02_assembly_%j.out

cd $SLURM_SUBMIT_DIR

# Run MEGAHIT assembly
eval "${ENV_CMD} bash snippets/02_assembly.sh"