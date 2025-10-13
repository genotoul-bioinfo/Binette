#!/bin/bash
#SBATCH --job-name=tutorial_05_binette
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --output=logs/tutorial_05_binette_%j.out
#SBATCH --error=logs/tutorial_05_binette_%j.err

# Load conda environment
source ~/miniconda3/bin/activate
conda activate binette_tutorial
cd $SLURM_SUBMIT_DIR

# Run Binette to refine and select best bins
bash snippets/05_binette.sh