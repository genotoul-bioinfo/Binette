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
eval "${ENV_CMD} bash snippets/01_download_dataset.sh"