[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/binette/README.html)  [![Anaconda-Server Badge](https://anaconda.org/bioconda/binette/badges/downloads.svg)](https://anaconda.org/bioconda/binette)

# Overview 

Binette is a fast and accurate binning refinement tool to constructs high quality MAGs from the output of multiple binning tools.

From the input bin sets, Binette constructs new hybrid bins. A bin can be seen as a set of contigs. When at least two bins overlap, meaning they share at least one contig, Binette utilizes basic set operations to create new bins.
- Intersection bin: This bin consists of the contigs that are shared by the overlapping bins. 
- Difference bin: This bin contains the contigs that are exclusively found in one bin and not present in the others.
- Union bin: The union bin includes all the contigs contained within the overlapping bins

It then uses checkm2 to assess bins quality to finally select the best bins possible.

Binette is inspired from the metaWRAP bin-refinement tool but it effectively solves all the problems from that very tool. 
- Enhanced Speed: Binette significantly improves the speed of the refinement process. It achieves this by launching the initial steps of checkm2, such as prodigal and diamond runs, only once on all contigs. These intermediate results are then utilized to assess the quality of any given bin, eliminating redundant computations and accelerating the refinement process.
- No Limit on Input Bin Sets: Unlike its predecessor, Binette is not constrained by the number of input bin sets. It can handle and process multiple bin sets simultaneously.
<!-- - Bin selection have been improved. It selects the best bins in a more accurate and elegant manner.
- It is easier to use. -->

# Installation

## With Bioconda

Binette can be esailly installed with conda 

```bash

conda install -c bioconda binette

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

Download checkm2 database

```

checkm2 database --download --path <checkm2/database/>
```

Finally install binette with pip

```
pip install .
```

Binette should be able to run :

```
binette -h
```


<!-- ## Running binette using a singularity image

Singularity version 3 or above must be installed. See [here](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps) how to install Singularity >=v3.

Git clone binette repository and build the singularity image. 

```
git clone https://github.com/genotoul-bioinfo/Binette
cd Binette
sudo singularity build binette.sif singularity_recipe
```

Then if the build succesfully finished, you should be able to run Binette:

```
singularity exec binette.sif binette -h
``` -->

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


# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker:

[binette issue tracker](https://github.com/genotoul-bioinfo/Binette/issues)

# Licence

This program is released as an open source software under the terms of [MIT License](https://forgemia.inra.fr/jean.mainguy/binette/-/raw/main/LICENSE).

