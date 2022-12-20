
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import sys
import logging
import os

import os 
from collections import Counter
from itertools import islice
import pandas as pd
import numpy as np
import logging
import logging
import os
import pyfastx

import itertools
import networkx as nx

# import pkg_resources
import pickle
from tqdm import tqdm

logging.basicConfig(level=logging.DEBUG,
                format='%(asctime)s %(levelname)s - %(message)s',
                datefmt="[%Y-%m-%d %H:%M:%S]")
                

def get_bin_graph(bins):
    logging.debug('start..')
    G = nx.Graph()
    logging.debug('add nodes..')
    G.add_nodes_from((b.id for b in bins))    

    logging.debug('Add edge')
    for i, (bin1, bin2) in tqdm(enumerate(itertools.combinations(bins, 2))):

            if bin1.overlaps_with(bin2):
                # logging.info(f"{bin1} overlaps with {bin2}")
                G.add_edge(bin1.id, bin2.id)
    return G


def select_best_bins(bins):
    logging.info(f'Building no bin graph from {len(bins)} ins')
    # G = get_bin_graph(bins)

    # nx.write_edgelist(G, "bin_graph_edglist")
    logging.info('SORTNG')
    sorted_bins = sorted(bins, key=lambda x: x.score, reverse=True)

    # remaining_bins = {b for b in bins}
    logging.info('SELECTING')
    selected_bins = []
    for b in sorted_bins:
        if b in bins:
            overlapping_bins = {b2 for b2 in bins if b.overlaps_with(b2)}
            bins -= overlapping_bins

            selected_bins.append(b)
            
            # neigbors_bins = nx.neighbors(G, b.id)
            # G.remove_nodes_from(list(neigbors_bins))

    logging.info(f'SELECTING {len(selected_bins)} bins')
    return selected_bins

def write_bin_info(bins, output):

    header = ['origin', "name", 'id', 'completeness', 'contamination', "score", 'size', 'contig_count',  'contigs']
    with open(output, "w") as fl:
        fl.write('\t'.join(header)+"\n")
        for bin_obj in bins:

            line = [bin_obj.origin, bin_obj.name, bin_obj.id,
                    bin_obj.completeness, bin_obj.contamination, bin_obj.score, 
                    bin_obj.length,
                     len(bin_obj.contigs), 
                     ";".join((str(c) for c in bin_obj.contigs))]

            fl.write("\t".join((str(e) for e in line)) + '\n')


with open("all_bin.p", "br") as fl:
    all_bins = pickle.load(fl)

logging.info(len(all_bins))

logging.info("SELECTION WITHOUT GRAPH ")

selected_bins = select_best_bins(all_bins)
logging.info("SELECTION WITHOUT GRAPH DONE")

logging.info('WRITING BINS')
write_bin_info(selected_bins, 'best_bins_no_graph.tsv')
