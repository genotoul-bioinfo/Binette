#!/bin/bash
#SBATCH --job-name=binette_tutorial_pipeline
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=00:10:00
#SBATCH --output=logs/master_tutorial_%j.out

# Master submission script for Binette tutorial pipeline
# This script submits all tutorial steps as dependent jobs

# Load environment
source ~/miniconda3/bin/activate

# Set working directory
cd $SLURM_SUBMIT_DIR

# Create logs directory if it doesn't exist
mkdir -p logs

# Print pipeline information
echo "=== Binette Tutorial Pipeline Master Script ==="
echo "Submitting complete tutorial pipeline as dependent jobs"
echo "Start time: $(date)"
echo ""

# Submit jobs with dependencies
echo "Submitting job chain..."

# Step 1: Download dataset
JOB1=$(sbatch --parsable sbatch_scripts/01_download_dataset.sbatch)
echo "Job 1 (Download): $JOB1"

# Step 2: Assembly (depends on step 1)
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 sbatch_scripts/02_assembly.sbatch)
echo "Job 2 (Assembly): $JOB2"

# Step 3: Read alignment (depends on step 2)
JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 sbatch_scripts/03_read_alignment.sbatch)
echo "Job 3 (Read alignment): $JOB3"

# Step 4: Binning (depends on step 3)
JOB4=$(sbatch --parsable --dependency=afterok:$JOB3 sbatch_scripts/04_binning.sbatch)
echo "Job 4 (Binning): $JOB4"

# Step 5: Binette (depends on step 4)
JOB5=$(sbatch --parsable --dependency=afterok:$JOB4 sbatch_scripts/05_binette.sbatch)
echo "Job 5 (Binette): $JOB5"

# Step 6: Analysis (depends on step 5)
JOB6=$(sbatch --parsable --dependency=afterok:$JOB5 sbatch_scripts/06_analysis.sbatch)
echo "Job 6 (Analysis): $JOB6"

echo ""
echo "=== Pipeline Submitted Successfully ==="
echo "Job dependency chain: $JOB1 -> $JOB2 -> $JOB3 -> $JOB4 -> $JOB5 -> $JOB6"
echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo "  watch 'squeue -u \$USER'"
echo ""
echo "Check job logs in logs/ directory:"
echo "  ls -la logs/"
echo ""
echo "Cancel entire pipeline if needed:"
echo "  scancel $JOB1 $JOB2 $JOB3 $JOB4 $JOB5 $JOB6"

echo "Master script completed at: $(date)"