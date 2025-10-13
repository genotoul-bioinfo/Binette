#!/bin/bash
#SBATCH --job-name=tutorial_02_assembly
#SBATCH --cpus-per-task=12
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --output=logs/tutorial_02_assembly_%j.out
#SBATCH --error=logs/tutorial_02_assembly_%j.err

# Load conda environment
source ~/miniconda3/bin/activate
conda activate binette_tutorial
cd $SLURM_SUBMIT_DIR

# Run MEGAHIT assembly
bash snippets/02_assembly.sh