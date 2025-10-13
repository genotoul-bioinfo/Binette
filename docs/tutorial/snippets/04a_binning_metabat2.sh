#!/bin/bash
# Generate depth file from BAM file for MetaBAT2 and MaxBin2
jgi_summarize_bam_contig_depths --outputDepth depth_Kickstart.txt Kickstart.bam

# Run MetaBAT2
mkdir -p metabat2_bins
metabat2 --inFile Kickstart.megahit/R1.contigs.fa --abdFile depth_Kickstart.txt --outFile metabat2_bins/metabat2 --numThreads 12 --seed 1