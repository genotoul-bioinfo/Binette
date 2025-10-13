#!/bin/bash
#SBATCH --job-name=tutorial_03_read_alignment
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --output=logs/tutorial_03_read_alignment_%j.out
#SBATCH --error=logs/tutorial_03_read_alignment_%j.err

# Load conda environment
source ~/miniconda3/bin/activate
conda activate binette_tutorial
cd $SLURM_SUBMIT_DIR

# Run BWA-MEM2 read alignment
bash snippets/03_read_alignment.sh