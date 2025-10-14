#!/bin/bash
#SBATCH --job-name=tutorial_06_analysis
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --time=00:30:00
#SBATCH --output=logs/tutorial_06_analysis_%j.out

cd $SLURM_SUBMIT_DIR

# Execute the analysis notebook and convert to HTML
eval "${ENV_CMD} jupyter nbconvert --to notebook --execute --inplace analyse_binette_result.ipynb"
