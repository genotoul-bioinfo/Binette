#!/bin/bash
# Download the Kickstart dataset using SRA toolkit
prefetch SRR5058924

# Convert SRA to paired FASTQ files with gzip compression
fastq-dump --defline-seq '@$ac_$sn/$ri' --defline-qual '+' --split-3 -O . --gzip SRR5058924/SRR5058924.sra

# Optional cleanup: remove the SRA file as it's no longer needed
rm -f SRR5058924/SRR5058924.sra