[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/binette/README.html)  [![Anaconda-Server Badge](https://anaconda.org/bioconda/binette/badges/downloads.svg)](https://anaconda.org/bioconda/binette)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/binette/badges/license.svg)](https://anaconda.org/bioconda/binette) 
[![Anaconda-Server Badge](https://anaconda.org/bioconda/binette/badges/version.svg)](https://anaconda.org/bioconda/binette)
[![PyPI version](https://badge.fury.io/py/Binette.svg)](https://badge.fury.io/py/Binette)
[![status](https://joss.theoj.org/papers/ad304709d59f1a51a31614393b09ba2b/status.svg)](https://joss.theoj.org/papers/ad304709d59f1a51a31614393b09ba2b)

[![Test Coverage](https://genotoul-bioinfo.github.io/Binette/coverage-badge.svg)](https://genotoul-bioinfo.github.io/Binette/) 
[![CI Status](https://github.com/genotoul-bioinfo/Binette/actions/workflows/binette_ci.yml/badge.svg)](https://github.com/genotoul-bioinfo/Binette/actions/workflows)
[![Documentation Status](https://readthedocs.org/projects/binette/badge/?version=latest)](https://binette.readthedocs.io/en/latest/?badge=latest)


# Binette 

Binette is a fast and accurate binning refinement tool to constructs high quality MAGs from the output of multiple binning tools.

From the input bin sets, Binette constructs new hybrid bins. A bin can be seen as a set of contigs. When at least two bins overlap, meaning they share at least one contig, Binette utilizes basic set operations to create new bins.
- Intersection bin: This bin consists of the contigs that are shared by the overlapping bins. 
- Difference bin: This bin contains the contigs that are exclusively found in one bin and not present in the others.
- Union bin: The union bin includes all the contigs contained within the overlapping bins

It then uses CheckM2 to assess bins quality to finally select the best bins possible.

Binette is inspired from the metaWRAP bin-refinement tool but it effectively solves all the problems from that very tool. 
- Enhanced Speed: Binette significantly improves the speed of the refinement process. It achieves this by launching the initial steps of CheckM2, such as Prodigal and Diamond runs, only once on all contigs. These intermediate results are then utilized to assess the quality of any given bin, eliminating redundant computations and accelerating the refinement process.
- No Limit on Input Bin Sets: Unlike its predecessor, Binette is not constrained by the number of input bin sets. It can handle and process multiple bin sets simultaneously.
<!-- - Bin selection have been improved. It selects the best bins in a more accurate and elegant manner.
- It is easier to use. -->

A comprehensive documentation of Binette is avalaible here: https://binette.readthedocs.io/

# Installation

## With Bioconda

Binette can be easilly installed with conda 

```bash
conda create -c bioconda -c defaults -c conda-forge -n binette binette
conda activate binette
```

Binette should be able to run :

```
binette -h
```


## From a conda environnement

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


## Downloading the CheckM2 database

Before using Binette, it is necessary to download the CheckM2 database:

```bash
checkm2 database --download --path <checkm2/database/>
```

Make sure to replace `<checkm2/database/>` with the desired path where you want to store the CheckM2 database.


# Usage 

## Input Formats

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

## Outputs

Binette results are stored in the `results` directory. You can specify a different directory using the `--outdir` option.

In this directory you will find:
- `final_bins_quality_reports.tsv`: This is a TSV (tab-separated values) file containing quality information about the final selected bins.
- `final_bins/`: This directory stores all the selected bins in fasta format.
- `temporary_files/`: This directory contains intermediate files. If you choose to use the `--resume` option, Binette will utilize files in this directory to prevent the recomputation of time-consuming steps.


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

# Help, feature requests and bug reporting

To report bugs, request new features, or seek help and support, please open an [issue](https://github.com/genotoul-bioinfo/Binette/issues). 


# Licence

This program is released as an open source software under the terms of [MIT License](LICENSE).

