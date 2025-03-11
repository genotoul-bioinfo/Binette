[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/binette/README.html)  [![Anaconda-Server Badge](https://anaconda.org/bioconda/binette/badges/downloads.svg)](https://anaconda.org/bioconda/binette)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/binette/badges/license.svg?cache-control=no-cache)](https://anaconda.org/bioconda/binette) 
[![Anaconda-Server Badge](https://anaconda.org/bioconda/binette/badges/version.svg?cache-control=no-cache)](https://anaconda.org/bioconda/binette)
[![PyPI version](https://badge.fury.io/py/Binette.svg?cache-control=no-cache)](https://badge.fury.io/py/Binette)
[![status](https://joss.theoj.org/papers/ad304709d59f1a51a31614393b09ba2b/status.svg)](https://joss.theoj.org/papers/ad304709d59f1a51a31614393b09ba2b)

[![Test Coverage](https://genotoul-bioinfo.github.io/Binette/coverage-badge.svg)](https://genotoul-bioinfo.github.io/Binette/) 
[![CI Status](https://github.com/genotoul-bioinfo/Binette/actions/workflows/binette_ci.yml/badge.svg)](https://github.com/genotoul-bioinfo/Binette/actions/workflows)
[![Documentation Status](https://readthedocs.org/projects/binette/badge/?version=latest)](https://binette.readthedocs.io/en/latest/?badge=latest)


# Binette

**Binette** is a fast and accurate binning refinement tool designed to construct high-quality MAGs from the output of multiple binning tools.  

### How It Works  

From the input bin sets, Binette constructs new **hybrid bins**. A **bin** is a set of contigs. When at least two bins overlap (i.e., share at least one contig), Binette applies fundamental set operations to generate new bins:  

- **Intersection Bin**: Contains contigs that are shared by overlapping bins.  
- **Difference Bin**: Includes contigs that are unique to one bin and not found in others.  
- **Union Bin**: Encompasses all contigs from the overlapping bins.  

Binette then evaluates bin quality using **CheckM2**, ensuring the selection of the best possible bins.  

### Why Binette?  

Binette is inspired by the **metaWRAP bin-refinement tool** but effectively addresses its limitations. Key improvements include:  

- **Enhanced Speed**  
  Binette significantly accelerates the refinement process by running **CheckM2's** initial steps (e.g., Prodigal and Diamond) only once for all contigs. These intermediate results are then reused for bin quality assessment, eliminating redundant computations.  

- **No Limit on Input Bin Sets**  
  Unlike metaWRAP, Binette supports **any number of input bin sets**, allowing seamless processing of multiple binning outputs.  

<!--  
- **Improved Bin Selection**  
  Binette selects the best bins in a more accurate and elegant manner.  

- **User-Friendly**  
  Designed for ease of use with a streamlined workflow.  
-->


> [!NOTE]
> For a detailed guide and tutorial, see the [Binette documentation](https://binette.readthedocs.io/). 


## Installation

### With Bioconda

Binette can be easilly installed with conda 

```bash
conda create -c bioconda -c defaults -c conda-forge -n binette binette
conda activate binette
```

Binette should be able to run :

```
binette -h
```


### From a conda environnement

Clone this repository: 
```
git clone https://github.com/genotoul-bioinfo/Binette
cd Binette
```

Then create a Conda environment using the `binette.yaml` file:
```
conda env create -n binette -f binette.yaml
conda activate binette 
```

Finally install Binette with pip

```
pip install .
```

Binette should be able to run :

```
binette -h
```


### Downloading the CheckM2 database

Before using Binette, it is necessary to download the CheckM2 database:

```bash
checkm2 database --download --path <checkm2/database/>
```

Make sure to replace `<checkm2/database/>` with the desired path where you want to store the CheckM2 database.


## Usage 

### Input Formats

Binette supports two input formats for bin sets: 

1. **Contig2bin Tables:** You can provide bin sets using contig2bin tables, which establish the relationship between each contig and its corresponding bin. In this format, you need to specify the `--contig2bin_tables` argument. 

For example, consider the following two `contig2bin_tables`:

- `bin_set1.tsv`:

    ```tsv
    contig_1   binA
    contig_8   binA
    contig_15  binB
    contig_9   binC
    ```
    
- `bin_set2.tsv`:

    ```tsv
    contig_1   bin.0
    contig_8   bin.0
    contig_15  bin.1
    contig_9   bin.2
    contig_10  bin.0
    ```
    
    The `binette` command to process this input would be:
    
    ```bash
    binette --contig2bin_tables bin_set1.tsv bin_set2.tsv --contigs assembly.fasta
    ```

2. **Bin Directories:** Alternatively, you can use bin directories, where each bin is represented by a separate FASTA file. For this format, you need to provide the `--bin_dirs` argument. Here's an example of two bin directories:

    ```
    bin_set1/
    ├── binA.fa: contains sequences of contig_1, contig_8
    ├── binB.fa: contains sequences of contig_15
    └── binC.fa: contains sequences of contig_9
    ```
    
    ```
    bin_set2/
    ├── binA.fa: contains sequences of contig_1, contig_8, contig_10
    ├── binB.fa: contains sequences of contig_15
    └── binC.fa: contains sequences of contig_9
    ```
    
    The `binette` command to process this input would be:
    
    ```bash
    binette --bin_dirs bin_set1 bin_set2 --contigs assembly.fasta
    ```

In both formats, the `--contigs` argument should specify a FASTA file containing all the contigs found in the bins. Typically, this file would be the assembly FASTA file used to generate the bins. In these exemple the `assembly.fasta` file should contain at least the five contigs mentioned in the `contig2bin_tables` files or in the bin fasta files: `contig_1`, `contig_8`, `contig_15`, `contig_9`, and `contig_10`.

### Outputs

Binette results are stored in the `results` directory. You can specify a different directory using the `--outdir` option.

In this directory you will find:
- **`final_bins_quality_reports.tsv`**: This is a TSV (tab-separated values) file containing quality information about the final selected bins.
- **`final_bins/`**: This directory stores all the selected bins in fasta format.
- **`input_bins_quality_reports/`**: A directory storing quality reports for the input bin sets, with files following the same structure as `final_bins_quality_reports.tsv`.
- **`temporary_files/`**: This directory contains intermediate files. If you choose to use the `--resume` option, Binette will utilize files in this directory to prevent the recomputation of time-consuming steps.


The `final_bins_quality_reports.tsv` file contains the following columns:
| Column Name         | Description                                                                                                  |
|---------------------|--------------------------------------------------------------------------------------------------------------|
| **bin_id**          | This column displays the unique ID of the bin.                                                             |
| **origin**          | Indicates the source or origin of the bin, specifying from which bin set it originates or the intermediate set operation that created it. |
| **name**            | The name of the bin.                                                                                        |
| **completeness**    | The completeness of the bin, determined by CheckM2.                                                         |
| **contamination**   | The contamination of the bin, determined by CheckM2.                                                       |
| **score**           | This column displays the computed score, which is calculated as: `completeness - contamination * weight`. You can customize the contamination weight using the `--contamination_weight` option. |
| **size**            | Represents the size of the bin in nucleotides.                                                              |
| **N50**             | Displays the N50 of the bin.                                                                                |
| **contig_count**    | The number of contigs contained within the bin.                                                             |

## Help, feature requests and bug reporting

To report bugs, request new features, or seek help and support, please open an [issue](https://github.com/genotoul-bioinfo/Binette/issues). 


## Licence

This tool is released as open source software under the terms of the [GNU General Public Licence](LICENSE).

## Citation  

Binette is a scientific software tool with a [published paper](https://joss.theoj.org/papers/10.21105/joss.06782) in the [Journal of Open Source Software](https://joss.theoj.org/). If you use Binette in academic research, please cite:  

> **Binette: a fast and accurate bin refinement tool to construct high-quality Metagenome Assembled Genomes.**  
> Mainguy et al., (2024).    
> *Journal of Open Source Software, 9(102), 6782.*    
> doi: [10.21105/joss.06782](https://doi.org/10.21105/joss.06782)   


Binette extensively uses **CheckM2**. If your work relies on Binette, consider citing:  

> **CheckM2: a rapid, scalable, and accurate tool for assessing microbial genome quality using machine learning.**   
> Chklovski, Alex, et al. (2023).  
> *Nature Methods, 20(8), 1203-1212.*    
> doi: [10.1038/s41592-023-01940-w](https://doi.org/10.1038/s41592-023-01940-w)


