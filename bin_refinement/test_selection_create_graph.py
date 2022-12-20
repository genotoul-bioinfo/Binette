
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


# def select_best_bins(bins):
#     logging.info(f'Building bin graph from {len(bins)} ins')
#     # G = get_bin_graph(bins)

#     # nx.write_edgelist(G, "bin_graph_edglist")
 
#     sorted_bins = sorted(bins, key=lambda x: x.score, reverse=True)

#     removed_bins = []
#     remaining_bins = []
#     selected_bins = []
#     for b in sorted_bins:
#         if b.id in remaining_bins:
#             overlapping_bins = [b2 for b2 in remaining_bins if b.overlaps_wit(b2)]

#             selected_bins.append(b)
            
#             # neigbors_bins = nx.neighbors(G, b.id)
#             # G.remove_nodes_from(list(neigbors_bins))
    
#     return selected_bins

logging.info("LOADIND PICKLE")

with open("all_bin.p", "br") as fl:
    all_bins = pickle.load(fl)

logging.info(len(all_bins))

logging.info("SORT")

sorted_bins = sorted(all_bins, key=lambda x: x.score, reverse=True)


logging.info("MAKE GRAPH")

G = get_bin_graph(all_bins)

logging.info("MAKE GRAPH DONE")
