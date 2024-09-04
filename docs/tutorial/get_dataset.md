## Obtaining Metagenomic Data for the Tutorial

### Using the ncezid-biome Datasets Tool

For this tutorial, we’ll use the "Kickstart" metagenome dataset from the [ncezid-biome datasets GitHub repository](https://github.com/ncezid-biome/). This dataset corresponds to sample [SAMN05024035](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRR5058924&o=acc_s%3Aa) and SRA [SRR5058924](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRR5058924&o=acc_s%3Aa).


We'll download the "Kickstart" dataset using the ncezid-biome datasets tool. You can find the tool and instructions on how to use it in their [GitHub repository](https://github.com/ncezid-biome/datasets?tab=readme-ov-file#edlb).

The tool called `uscdc-datasets-sars-cov-2` on bioconda is part of the Conda environment created in the [previous section](./set_environment.md). 


#### Download the Kickstart Dataset

Once the tool is installed, you can download the "Kickstart" dataset with the following steps:

1. **Download the coal-metagenomics table** from the GitHub repository:
   
   ```{code-block} bash
   wget https://raw.githubusercontent.com/ncezid-biome/datasets/master/datasets/coal-metagenomics.tsv
   ```

2. **Select the relevant line** corresponding to the "Kickstart" dataset (SRR5058924) by extracting the header and the specific entry:

   ```{code-block} bash
   # Select the header of the table
   head -n7 coal-metagenomics.tsv > coal-metagenomics_Kickstart_only.tsv
   
   # Append the relevant line for the Kickstart dataset
   grep SRR5058924 coal-metagenomics.tsv >> coal-metagenomics_Kickstart_only.tsv
   ```

3. **Run the dataset download** using the `GenFSGopher.pl` script:

   ```{code-block} bash
   GenFSGopher.pl --numcpus 12 --compressed --outdir coal-metagenomics coal-metagenomics_Kickstart_only.tsv
   ```


   :::{admonition} ⌛ Expected Time
   :class: note

   This process takes approximately 16 minutes to complete.
   :::

#### Directory Structure

After downloading, your directory structure should look like this:

```{code-block} text
├── coal-metagenomics_Kickstart_only.tsv
└── data
    ├── in.tsv
    ├── Kickstart_1.fastq.gz
    ├── Kickstart_1.fastq.sha256
    ├── Kickstart_2.fastq.gz
    ├── Kickstart_2.fastq.sha256
    ├── Makefile
    ├── prefetch.done
    ├── sha256sum.log
    ├── SRR5058924
    │   └── SRR5058924.sra
    └── tree.dnd
```

In the next section, will assemble the two reads files to obtain an assembly of the dataset:
- `data/Kickstart_1.fastq.gz`
- `data/Kickstart_2.fastq.gz`


:::{admonition} 🧹 Cleaning Tip
:class: tip

You can remove the SRA file `data/SRR5058924/SRR5058924.sra` as it is no longer needed; we will use only the FASTQ files. To remove it, run:

```{code-block} bash
rm data/SRR5058924/SRR5058924.sra
:::

```{note}
Alternatively, you can download the data using the SRA Toolkit, which is what the ncezid-biome tool uses in the background. 
Note that ncezid-biome tool provides additional checksum verification to ensure data integrity.
You can retrieve the data with the following commands after installing the SRA Toolkit (e.g., via Conda: [sra-tools on Anaconda](https://anaconda.org/bioconda/sra-tools)):
```{code-block} bash
prefetch SRR5058924
fastq-dump --defline-seq '@$ac_$sn/$ri' --defline-qual '+' --split-3 -O . SRR5058924.sra
```

