#!/bin/bash
#SBATCH --job-name=tutorial_05_binette
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --output=logs/tutorial_05_binette_%j.out
#SBATCH --error=logs/tutorial_05_binette_%j.err

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
required_dirs=("metabat2_bins" "maxbin2_bins" "concoct_bins" "semibin2_output")
for dir in "${required_dirs[@]}"; do
    if [ ! -d "$dir" ]; then
        echo "ERROR: Required directory $dir not found. Please run 04_binning.sbatch first."
        exit 1
    fi
done

# Execute the Binette script
echo "=== Running Binette Bin Refinement ==="
bash snippets/05_binette.sh

# Check results
if [ $? -eq 0 ] && [ -d "binette_output" ]; then
    echo ""
    echo "=== Binette completed successfully ==="
    echo "Final bin statistics:"
    echo "Number of final bins: $(ls binette_output/final_bin_set/*.fa 2>/dev/null | wc -l)"
    echo "Binette output directory:"
    ls -la binette_output/
    
    echo ""
    echo "Results summary from bin_info.tsv:"
    if [ -f "binette_output/bin_info.tsv" ]; then
        echo "Columns: $(head -1 binette_output/bin_info.tsv)"
        echo "Number of bins: $(tail -n +2 binette_output/bin_info.tsv | wc -l)"
        echo ""
        echo "Quality summary:"
        tail -n +2 binette_output/bin_info.tsv | cut -f3,4 | sort -k1,1nr | head -5
    fi
    
    echo ""
    echo "Total disk usage:"
    du -sh binette_output/
else
    echo ""
    echo "=== Binette failed ==="
    exit 1
fi

echo "Job completed at: $(date)"