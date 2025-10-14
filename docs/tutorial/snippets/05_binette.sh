#!/bin/bash
# Run Binette to refine and improve bins
binette --bin_dirs maxbin2_bins/ metabat2_bins/ semibin2_output/output_bins/ concoct_bins/ \
        -c Kickstart.megahit/R1.contigs.fa \
         -t 12 -o binette_output