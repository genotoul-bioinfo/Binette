## Assemble the Reads

We will use **MEGAHIT** to assemble the reads from our dataset. Run the following command:

```{code-block} bash
megahit -1 coal-metagenomics/Kickstart_1.fastq.gz \
        -2 coal-metagenomics/Kickstart_2.fastq.gz \
        --out-dir Kickstart.megahit --out-prefix R1 --num-cpu-threads 12
```

:::{admonition} ⌛ Expected Time
:class: note

This process takes approximately 28 minutes to complete.
:::


```{admonition} Note
:class: note

You can also use **SPAdes** for assembly. It generally performs better than MEGAHIT but takes longer and requires more memory. Refer to the CAMI benchmark for a detailed comparison.
```

## Align the Reads Over the Assembly

To get coverage information, we first need to map the reads back to the assembly.

```{code-block} bash
# Create a directory for the alignments
mkdir -p alignments_bwa/

# Index the contigs file using BWA-MEM2
bwa-mem2 index Kickstart.megahit/R1.contigs.fa -p Kickstart.megahit/R1.contigs.fa

# Map reads back to the assembly, convert to BAM format, and sort
bwa-mem2 mem -t 12 Kickstart.megahit/R1.contigs.fa coal-metagenomics/Kickstart_*.fastq.gz | \
samtools view -@ 12 -bS - | \
samtools sort -@ 12 - -o alignments_bwa/Kickstart.bam

# Index the BAM file
samtools index alignments_bwa/Kickstart.bam
```


:::{admonition} ⌛ Expected Time
:class: note

This process takes approximately 12 minutes to complete.
:::
