
## Set tutorial environment

We will download necessary tool in a dedicated conda envrionnement.


<!-- #region -->
Let's create a directory to run the tutorial:


```bash

mamba env create -f binette_tutorial_env.yaml -n binette_tuto

```
<!-- #endregion -->

<!-- #region -->
## Get the Data

### Using ncezid-biome datasets tool

I downloaded the metagenome Kickstart from the above dataset (SAMN05024035) that correspond to this sra SRR5058924 https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRR5058924&o=acc_s%3Aa


We will donwload the data of the Kickstart (SAMN05024035) dataset this repository that https://github.com/ncezid-biome/datasets?tab=readme-ov-file#edlb

We had use conda as detailed here https://github.com/ncezid-biome/datasets/blob/master/INSTALL.md#conda 

Now we can download the Kickstart dataset with the folowing commands. 

We first download the coal-metagenomic table from the  github repository : https://github.com/ncezid-biome/datasets/blob/master/datasets/coal-metagenomics.tsv
ANd just select the line corresponding to the Kickstart dataset.


<!-- #endregion -->

```bash
# download the coal-metagenomic tsv file from the github repository
wget https://raw.githubusercontent.com/ncezid-biome/datasets/master/datasets/coal-metagenomics.tsv

# select the header of the table as it is necessary for the download

head -n7 coal-metagenomics.tsv > coal-metagenomics_Kickstart_only.tsv
grep SRR5058924  coal-metagenomics.tsv >> coal-metagenomics_Kickstart_only.tsv

GenFSGopher.pl --numcpus 12 --compressed --outdir  coal-metagenomics coal-metagenomics.tsv

```

It takes around 16min to run

You should hae the folowing structure
```
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

```{tip}
You can remove the SRA file `data/SRR5058924/SRR5058924.sra` as we do not need it anymore as we will exclusively use the fastq files. with `rm data/SRR5058924/SRR5058924.sra`
```

```{note}
You can also download the data using SRA toolkit which what the tool does in the background but add some check sum to ensure data integrity. After instaling sra toolkit (with conda for example : https://anaconda.org/bioconda/sra-tools) you can run the two commands folowing commands to retrived the data: `prefetch  SRR5058924` and `fastq-dump --defline-seq '@$ac_$sn/$ri' --defline-qual '+' --split-3 -O . SRR5058924.sra` 
```