#!/bin/bash
# Download the coal-metagenomics table from GitHub
wget https://raw.githubusercontent.com/ncezid-biome/datasets/master/datasets/coal-metagenomics.tsv

# Select the header of the table
head -n8 coal-metagenomics.tsv > coal-metagenomics_Kickstart_only.tsv

# Append the relevant line for the Kickstart dataset
grep SRR5058924 coal-metagenomics.tsv >> coal-metagenomics_Kickstart_only.tsv

# Run the dataset download using the GenFSGopher.pl script
GenFSGopher.pl --numcpus 12 --compressed --outdir coal-metagenomics coal-metagenomics_Kickstart_only.tsv

# Optional cleanup: remove the SRA file as it's no longer needed
rm -f coal-metagenomics/SRR5058924/SRR5058924.sra