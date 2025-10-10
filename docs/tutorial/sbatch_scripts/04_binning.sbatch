#!/bin/bash
#SBATCH --job-name=tutorial_04_binning
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --time=04:00:00
#SBATCH --output=logs/tutorial_04_binning_%j.out
#SBATCH --error=logs/tutorial_04_binning_%j.err

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
if [ ! -f "Kickstart.bam" ]; then
    echo "ERROR: BAM file not found. Please run 03_read_alignment.sbatch first."
    exit 1
fi

# Execute all binning scripts in sequence
echo "=== Running MetaBAT2 Binning ==="
bash snippets/04a_binning_metabat2.sh

if [ $? -ne 0 ]; then
    echo "ERROR: MetaBAT2 binning failed"
    exit 1
fi

echo ""
echo "=== Running MaxBin2 Binning ==="
bash snippets/04b_binning_maxbin2.sh

if [ $? -ne 0 ]; then
    echo "ERROR: MaxBin2 binning failed"
    exit 1
fi

echo ""
echo "=== Running CONCOCT Binning ==="
bash snippets/04c_binning_concoct.sh

if [ $? -ne 0 ]; then
    echo "ERROR: CONCOCT binning failed"
    exit 1
fi

echo ""
echo "=== Running SemiBin2 Binning ==="
bash snippets/04d_binning_semibin2.sh

if [ $? -ne 0 ]; then
    echo "ERROR: SemiBin2 binning failed"
    exit 1
fi

# Check results
echo ""
echo "=== Binning completed successfully ==="
echo "Binning results summary:"
echo "MetaBAT2 bins: $(ls metabat2_bins/*.fa 2>/dev/null | wc -l)"
echo "MaxBin2 bins: $(ls maxbin2_bins/*.fasta 2>/dev/null | wc -l)"
echo "CONCOCT bins: $(ls concoct_bins/*.fa 2>/dev/null | wc -l)"
echo "SemiBin2 bins: $(ls semibin2_output/output_bins/*.fa 2>/dev/null | wc -l)"

echo ""
echo "Disk usage by binning results:"
du -sh metabat2_bins/ maxbin2_bins/ concoct_bins/ semibin2_output/ 2>/dev/null || echo "Some directories not found"

echo "Job completed at: $(date)"