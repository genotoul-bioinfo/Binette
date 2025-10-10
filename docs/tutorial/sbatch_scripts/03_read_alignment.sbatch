#!/bin/bash
#SBATCH --job-name=tutorial_03_read_alignment
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --output=logs/tutorial_03_read_alignment_%j.out
#SBATCH --error=logs/tutorial_03_read_alignment_%j.err

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
if [ ! -f "Kickstart.megahit/R1.contigs.fa" ]; then
    echo "ERROR: Assembly not found. Please run 02_assembly.sbatch first."
    exit 1
fi

# Execute the read alignment script
echo "=== Running BWA Read Alignment ==="
bash snippets/03_read_alignment.sh

# Check results
if [ $? -eq 0 ] && [ -f "Kickstart.bam" ]; then
    echo ""
    echo "=== Read alignment completed successfully ==="
    echo "BAM file statistics:"
    samtools flagstat Kickstart.bam
    echo ""
    echo "BAM file size: $(du -h Kickstart.bam)"
    echo "Index file: $(ls -la Kickstart.bam.bai 2>/dev/null || echo "No index found")"
else
    echo ""
    echo "=== Read alignment failed ==="
    exit 1
fi

echo "Job completed at: $(date)"