#!/bin/bash
# Assemble the reads using MEGAHIT
megahit -1 coal-metagenomics/Kickstart_1.fastq.gz \
        -2 coal-metagenomics/Kickstart_2.fastq.gz \
        --out-dir Kickstart.megahit --out-prefix R1 --num-cpu-threads 12