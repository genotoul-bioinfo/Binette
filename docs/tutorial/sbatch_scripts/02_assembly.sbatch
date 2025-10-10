#!/bin/bash
#SBATCH --job-name=tutorial_02_assembly
#SBATCH --cpus-per-task=12
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --output=logs/tutorial_02_assembly_%j.out
#SBATCH --error=logs/tutorial_02_assembly_%j.err

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

# Check prerequisites
if [ ! -f "coal-metagenomics/Kickstart_1.fastq.gz" ]; then
    echo "ERROR: Input files not found. Please run 01_download_dataset.sbatch first."
    exit 1
fi

# Execute the assembly script
echo "=== Running MEGAHIT Assembly ==="
bash snippets/02_assembly.sh

# Check results
if [ $? -eq 0 ] && [ -f "Kickstart.megahit/R1.contigs.fa" ]; then
    echo ""
    echo "=== Assembly completed successfully ==="
    echo "Assembly statistics:"
    echo "Number of contigs: $(grep -c '^>' Kickstart.megahit/R1.contigs.fa)"
    echo "Assembly size: $(du -h Kickstart.megahit/R1.contigs.fa)"
    echo "Assembly directory:"
    ls -la Kickstart.megahit/
else
    echo ""
    echo "=== Assembly failed ==="
    exit 1
fi

echo "Job completed at: $(date)"