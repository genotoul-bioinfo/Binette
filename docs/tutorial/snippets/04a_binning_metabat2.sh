#!/bin/bash
# Generate depth file from BAM file for MetaBAT2 and MaxBin2
jgi_summarize_bam_contig_depths --outputDepth depth_Kickstart.txt alignments_bwa/Kickstart.bam

# Run MetaBAT2
metabat2 --inFile Kickstart.megahit/R1.contigs.fa --abdFile depth_Kickstart.txt --outFile metabat2/metabat2 --numThreads 12 --seed 1