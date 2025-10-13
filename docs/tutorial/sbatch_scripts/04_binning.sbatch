#!/bin/bash
#SBATCH --job-name=tutorial_04_binning
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --time=04:00:00
#SBATCH --output=logs/tutorial_04_binning_%j.out
#SBATCH --error=logs/tutorial_04_binning_%j.err

cd $SLURM_SUBMIT_DIR

# Run all four binning tools sequentially
${ENV_CMD:-source ~/miniconda3/bin/activate && conda activate binette_tutorial &&} bash snippets/04a_binning_metabat2.sh
${ENV_CMD:-source ~/miniconda3/bin/activate && conda activate binette_tutorial &&} bash snippets/04b_binning_maxbin2.sh
${ENV_CMD:-source ~/miniconda3/bin/activate && conda activate binette_tutorial &&} bash snippets/04c_binning_concoct.sh
${ENV_CMD:-source ~/miniconda3/bin/activate && conda activate binette_tutorial &&} bash snippets/04d_binning_semibin2.sh