#!/bin/bash
# Run MaxBin2 using the depth file from MetaBAT2
mkdir -p maxbin2_bins
run_MaxBin.pl -contig Kickstart.megahit/R1.contigs.fa \
                -abund depth_Kickstart.txt -thread 12 -out maxbin2_bins/maxbin2