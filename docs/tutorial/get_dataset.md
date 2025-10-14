## Obtaining Metagenomic Data for the Tutorial

### Using the SRA Toolkit

For this tutorial, we'll use the "Kickstart" metagenome dataset which corresponds to sample [SAMN05024035](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRR5058924&o=acc_s%3Aa) and SRA [SRR5058924](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRR5058924&o=acc_s%3Aa).

We'll download the dataset directly using the SRA toolkit. The SRA tools (`sra-tools` package) are included in the Conda environment created in the [previous section](./set_environment.md).

#### Download the Kickstart Dataset

You can download the "Kickstart" dataset with the following commands:

```{include} snippets/01_download_dataset.sh
:code: bash
```

:::{admonition} SRA Toolkit Information
:class: tip

The SRA toolkit provides direct access to sequencing data from NCBI's Sequence Read Archive. The `prefetch` command downloads the SRA file locally, and `fastq-dump` converts it to standard FASTQ format with proper paired-end splitting and gzip compression for efficiency.


You can remove the SRA file `SRR5058924/SRR5058924.sra` as it is no longer needed after conversion to FASTQ files. To remove it run:

```{code-block} bash
rm SRR5058924/SRR5058924.sra
```

:::


:::{admonition} ⌛ Expected Time
:class: note

This process takes approximately 5-10 minutes to complete.
:::

#### Directory Structure

After downloading, your directory structure should look like this:

```{code-block} text
├── SRR5058924/
├── SRR5058924_1.fastq.gz
└── SRR5058924_2.fastq.gz
```

The `prefetch` command downloads the SRA file to the `SRR5058924/` directory, and `fastq-dump` converts it to paired FASTQ files with gzip compression. The SRA file is automatically cleaned up after conversion.

In the next section, we will assemble the two reads files to obtain an assembly of the dataset:
- `SRR5058924_1.fastq.gz` (forward reads)  
- `SRR5058924_2.fastq.gz` (reverse reads)


