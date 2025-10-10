# SBATCH Scripts for Binette Tutorial

This directory contains SLURM SBATCH scripts to run the Binette tutorial on HPC clusters. These scripts use the same code snippets as the documentation, ensuring perfect synchronization between tutorial documentation and cluster execution.

## Quick Start

1. **Edit SBATCH headers** in each script to match your cluster requirements (partition, QoS, etc.)
2. **Submit the entire pipeline**:
   ```bash
   cd /path/to/binette/docs/tutorial
   sbatch sbatch_scripts/submit_tutorial_pipeline.sbatch
   ```

## Individual Script Execution

If you prefer to run steps individually:

```bash
# Step 1: Download dataset (prerequisite for all other steps)
sbatch sbatch_scripts/01_download_dataset.sbatch

# Step 2: Assembly (requires dataset)
sbatch sbatch_scripts/02_assembly.sbatch

# Step 3: Read alignment (requires assembly)
sbatch sbatch_scripts/03_read_alignment.sbatch

# Step 4: Binning (requires alignment)
sbatch sbatch_scripts/04_binning.sbatch

# Step 5: Binette (requires binning results)
sbatch sbatch_scripts/05_binette.sbatch
```

## Script Details

### Resource Requirements

| Script | CPUs | Memory | Time | Notes |
|--------|------|--------|------|-------|
| `01_download_dataset.sbatch` | 12 | 16GB | 30min | Network intensive |
| `02_assembly.sbatch` | 12 | 32GB | 2h | Memory intensive |
| `03_read_alignment.sbatch` | 8 | 16GB | 1h | CPU intensive |
| `04_binning.sbatch` | 16 | 64GB | 4h | Most resource intensive |
| `05_binette.sbatch` | 8 | 32GB | 2h | Final refinement |

### Environment Setup

All scripts assume:
- Miniconda3 installed at `~/miniconda3/`
- Conda environment named `binette_tutorial` exists
- Environment created from `binette_tutorial_env.yaml`

To create the environment:
```bash
conda env create -f binette_tutorial_env.yaml
```

### Directory Structure

Scripts create the following structure:
```
tutorial_workdir/
├── logs/                           # SLURM job logs
├── coal-metagenomics/             # Downloaded dataset
├── Kickstart.megahit/            # Assembly results
├── Kickstart.bam                 # Read alignment
├── metabat2_bins/                 # MetaBAT2 bins
├── maxbin2_bins/                  # MaxBin2 bins
├── concoct_bins/                  # CONCOCT bins
├── semibin2_output/               # SemiBin2 bins
└── binette_output/                # Final Binette results
```

## Code Synchronization

These SBATCH scripts source the same shell snippets used in the Sphinx documentation:

- `snippets/01_download_dataset.sh`
- `snippets/02_assembly.sh`
- `snippets/03_read_alignment.sh`
- `snippets/04a_binning_metabat2.sh`
- `snippets/04b_binning_maxbin2.sh`
- `snippets/04c_binning_concoct.sh`
- `snippets/04d_binning_semibin2.sh`
- `snippets/05_binette.sh`

This ensures that:
- ✅ Documentation and cluster scripts are always synchronized
- ✅ Changes to tutorial commands automatically propagate to cluster scripts
- ✅ No code duplication between markdown and SBATCH files
- ✅ Single source of truth for tutorial commands

## Monitoring and Troubleshooting

### Check Job Status
```bash
# Monitor all your jobs
squeue -u $USER

# Watch job progress
watch 'squeue -u $USER'

# Check specific job details
scontrol show job <JOBID>
```

### Check Logs
```bash
# View recent job logs
ls -la logs/

# Check specific step logs
tail -f logs/tutorial_02_assembly_<JOBID>.out
tail -f logs/tutorial_02_assembly_<JOBID>.err
```

### Common Issues

1. **Environment not found**: Ensure conda environment exists
   ```bash
   conda env list
   conda env create -f binette_tutorial_env.yaml
   ```

2. **Permission errors**: Check file permissions and paths
3. **Resource limits**: Adjust SBATCH headers for your cluster
4. **Missing dependencies**: Ensure all tools are installed in conda environment

### Cancel Jobs
```bash
# Cancel specific job
scancel <JOBID>

# Cancel all your jobs
scancel -u $USER

# Cancel entire pipeline (if using master script)
# Job IDs will be printed by submit_tutorial_pipeline.sbatch
scancel <JOB1> <JOB2> <JOB3> <JOB4> <JOB5>
```

## Customization

### Cluster-Specific Settings

Edit SBATCH headers in each script:
```bash
#SBATCH --partition=<your_partition>
#SBATCH --qos=<your_qos>
#SBATCH --account=<your_account>
```

### Resource Adjustments

Modify CPU and memory allocations based on your data size:
```bash
#SBATCH --cpus-per-task=<cores>
#SBATCH --mem=<memory>
#SBATCH --time=<time_limit>
```

### Environment Path

Update conda activation path if different:
```bash
source /path/to/your/miniconda3/bin/activate
```

## Testing

To test scripts without running full pipeline:

1. **Dry run**: Check SBATCH syntax
   ```bash
   sbatch --test-only sbatch_scripts/01_download_dataset.sbatch
   ```

2. **Small dataset**: Use reduced parameters for testing
3. **Interactive session**: Test commands interactively first
   ```bash
   srun --pty bash
   source ~/miniconda3/bin/activate
   conda activate binette_tutorial
   bash snippets/01_download_dataset.sh
   ```

## Contributing

When updating tutorial commands:

1. ✅ Update the snippet files in `snippets/`
2. ✅ SBATCH scripts automatically use updated snippets
3. ✅ Documentation automatically includes updated snippets
4. ✅ Test on cluster to ensure compatibility

This ensures perfect synchronization across all tutorial formats!