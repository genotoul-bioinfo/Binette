#!/bin/bash
# Run SemiBin2 with single_easy_bin command
SemiBin2 single_easy_bin -i Kickstart.megahit/R1.contigs.fa \
                            -b alignments_bwa/Kickstart.bam \
                            -o semibin2/ -p 12