# SBATCH Scripts for Binette Tutorial

SLURM scripts to run the Binette tutorial on HPC clusters. Scripts use the same code snippets as the documentation for perfect synchronization.

## Usage

### With Conda (default)
Create conda environment:
```bash
conda env create -f ../binette_tutorial_env.yaml
```

Submit pipeline (uses conda by default):
```bash
sbatch submit_tutorial_pipeline.sbatch
```

### With Apptainer
Set the ENV_CMD variable to use Apptainer:
```bash
export ENV_CMD="apptainer exec /path/to/binette.sif"
sbatch submit_tutorial_pipeline.sbatch
```

### Individual steps
```bash
sbatch 01_download_dataset.sbatch
sbatch 02_assembly.sbatch  
sbatch 03_read_alignment.sbatch
sbatch 04_binning.sbatch
sbatch 05_binette.sbatch
```

## Configuration

### Environment Command
Set `ENV_CMD` environment variable to control execution environment:

- **Conda** (default): `ENV_CMD="source ~/miniconda3/bin/activate && conda activate binette_tutorial &&"`
- **Apptainer**: `ENV_CMD="apptainer exec /path/to/image.sif"`
- **Custom**: Any command prefix you need

### Examples
```bash
# Use different conda path
export ENV_CMD="source /opt/conda/bin/activate && conda activate binette_tutorial &&"

# Use Apptainer with bind mounts
export ENV_CMD="apptainer exec --bind /data:/data /path/to/binette.sif"

# Use module system + conda
export ENV_CMD="module load conda && conda activate binette_tutorial &&"
```

## Notes

- Edit SBATCH headers to match your cluster (partition, account, etc.)
- All error checking is handled by the snippet scripts themselves