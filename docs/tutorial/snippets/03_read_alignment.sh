#!/bin/bash
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