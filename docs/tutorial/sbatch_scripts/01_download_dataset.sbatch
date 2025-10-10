#!/bin/bash
#SBATCH --job-name=tutorial_01_download
#SBATCH --cpus-per-task=12
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --output=logs/tutorial_01_download_%j.out
#SBATCH --error=logs/tutorial_01_download_%j.err

# Load environment
source ~/miniconda3/bin/activate
conda activate binette_tutorial

# Set working directory
cd $SLURM_SUBMIT_DIR

# Print job information
echo "Job ID: $SLURM_JOB_ID"
echo "Job name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

# Execute the download script
echo "=== Running Dataset Download ==="
bash snippets/01_download_dataset.sh

# Check results
if [ $? -eq 0 ]; then
    echo ""
    echo "=== Download completed successfully ==="
    echo "Files created:"
    ls -la coal-metagenomics/
    echo ""
    echo "Dataset info:"
    wc -l coal-metagenomics_Kickstart_only.tsv
else
    echo ""
    echo "=== Download failed ==="
    exit 1
fi

echo "Job completed at: $(date)"