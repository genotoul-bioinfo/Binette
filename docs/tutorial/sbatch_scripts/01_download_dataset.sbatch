#!/bin/bash
#SBATCH --job-name=tutorial_01_download
#SBATCH --cpus-per-task=12
#SBATCH --mem=16GB
#SBATCH --time=00:30:00
#SBATCH --output=logs/tutorial_01_download_%j.out

# Set working directory and create logs
cd $SLURM_SUBMIT_DIR
mkdir -p logs

# Run the download script
${ENV_CMD:-source ~/miniconda3/bin/activate && conda activate binette_tutorial &&} bash snippets/01_download_dataset.sh