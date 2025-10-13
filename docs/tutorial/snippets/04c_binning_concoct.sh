#!/bin/bash
# Run CONCOCT binning

# Create directory
mkdir -p concoct/

# Cut up the FASTA file into chunks for processing
cut_up_fasta.py Kickstart.megahit/R1.contigs.fa --chunk_size 10000 \
                --overlap_size 0 --merge_last \
                --bedfile concoct/contigs_10K.bed > concoct/contigs_10K.fa

# Generate the coverage table from the BAM file
concoct_coverage_table.py concoct/contigs_10K.bed Kickstart.bam > concoct/coverage_table.tsv

# Run CONCOCT with the composition and coverage files
concoct --composition_file concoct/contigs_10K.fa \
        --coverage_file concoct/coverage_table.tsv \
        --basename concoct/bins --threads 12

# Merge the clustering results and extract bins
merge_cutup_clustering.py concoct/bins_clustering_gt1000.csv > concoct/clustering_merge.csv

mkdir -p concoct_bins

extract_fasta_bins.py Kickstart.megahit/R1.contigs.fa concoct/clustering_merge.csv --output_path concoct_bins