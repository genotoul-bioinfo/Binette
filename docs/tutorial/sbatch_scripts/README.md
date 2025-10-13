# SBATCH Scripts for Binette Tutorial

SLURM scripts to run the Binette tutorial on HPC clusters. Scripts use the same code snippets as the documentation for perfect synchronization.

## Usage

Create conda environment:
```bash
conda env create -f ../binette_tutorial_env.yaml
```

Submit entire pipeline:
```bash
sbatch submit_tutorial_pipeline.sbatch
```

Or submit individual steps:
```bash
sbatch 01_download_dataset.sbatch
sbatch 02_assembly.sbatch  
sbatch 03_read_alignment.sbatch
sbatch 04_binning.sbatch
sbatch 05_binette.sbatch
```

## Notes

- Edit SBATCH headers to match your cluster (partition, account, etc.)
- Scripts assume `~/miniconda3/` conda installation
- All error checking is handled by the snippet scripts themselves