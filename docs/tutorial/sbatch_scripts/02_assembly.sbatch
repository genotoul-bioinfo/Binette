#!/bin/bash
#SBATCH --job-name=tutorial_02_assembly
#SBATCH --cpus-per-task=12
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --output=logs/tutorial_02_assembly_%j.out
#SBATCH --error=logs/tutorial_02_assembly_%j.err

cd $SLURM_SUBMIT_DIR

# Run MEGAHIT assembly
${ENV_CMD:-source ~/miniconda3/bin/activate && conda activate binette_tutorial &&} bash snippets/02_assembly.sh