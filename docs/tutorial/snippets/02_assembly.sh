#!/bin/bash
# Assemble the reads using MEGAHIT
megahit -1 SRR5058924_1.fastq.gz \
        -2 SRR5058924_2.fastq.gz \
        --out-dir Kickstart.megahit --out-prefix R1 --num-cpu-threads 12