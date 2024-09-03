

## Assemble the reads

We will use megahit to assemble the reads

```bash

cd /home/jmainguy/Analysis/Binette/Binette_tutorial/ncezid-biome_datasets/exec_tutorial_jupyter
```

```bash

megahit -1 coal-metagenomics/Kickstart_1.fastq.gz -2 coal-metagenomics/Kickstart_2.fastq.gz --out-dir Kickstart.megahit --out-prefix R1 --num-cpu-threads 12

```


This take 27m49,879s 

```{note}
We can use spade as well. It performs generally better that megahit but is generally longer and consume more memory than megahit. See cami benchmark ???  
```



## Align the reads over the assembly 

First we need to map the reads back against the assembly to get coverage information

```bash

mkdir -p alignments_bwa/

bwa-mem2 index Kickstart.megahit/R1.contigs.fa -p Kickstart.megahit/R1.contigs.fa

bwa-mem2 mem -t 12  Kickstart.megahit/R1.contigs.fa  coal-metagenomics/Kickstart_*.fastq.gz | samtools view -@ 12 -bS - | samtools sort -@ 12  - -o alignments_bwa/Kickstart.bam

samtools index alignments_bwa/Kickstart.bam 

```
<!-- #region -->
This take around   12 minutes