"""


"""

import pyfastx 
import logging

def parse_fasta_file(fasta_file):
    fa = pyfastx.Fasta(fasta_file, build_index=True)
    return fa


def make_contig_index(contigs):
    contig_to_index = {contig:index for index, contig in enumerate(contigs)}
    index_to_contig = {index:contig for contig, index  in contig_to_index.items() }
    return contig_to_index, index_to_contig


def apply_contig_index(contig_to_index, contig_to_info):
    return {contig_to_index[contig]:info for contig, info in contig_to_info.items()}

