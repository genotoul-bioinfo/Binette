#!/bin/bash
#SBATCH --job-name=tutorial_04_binning
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --output=logs/tutorial_04_binning_%j.out

cd $SLURM_SUBMIT_DIR

# Run all four binning tools sequentially
# Run binning tools
eval "${ENV_CMD} bash snippets/04a_binning_metabat2.sh"
eval "${ENV_CMD} bash snippets/04b_binning_maxbin2.sh"
eval "${ENV_CMD} bash snippets/04c_binning_concoct.sh"
eval "${ENV_CMD} bash snippets/04d_binning_semibin2.sh"
