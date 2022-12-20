
import sys
import logging
import bin_manager
import file_manager
import os
import cds
import diamond
# import bin_quality


logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt="[%Y-%m-%d %H:%M:%S]")

faa_file = 'tmp_head.faa'
t=30

logging.info('Parse faa file.')
contig_to_genes = cds.parse_faa_file(faa_file)

# from genes extract contig cds metadata
logging.info('Compute cds metadata.')

contig_to_genes_ss = {}
for i, (c, g) in enumerate(contig_to_genes.items()):
    if i == 10000:
        break
    contig_to_genes_ss[c] = g

contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length = cds.get_contig_cds_metadata(contig_to_genes_ss, t) 
