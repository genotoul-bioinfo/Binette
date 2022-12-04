"""


"""

import logging
import os
import pyfastx

import itertools
import networkx as nx

class Bin:
    counter = 0
    def __init__(self, contigs, origin, name):
        Bin.counter += 1 

        self.origin = origin
        self.name = name
        self.id = Bin.counter
        self.contigs = set(contigs)
        self.length = None

    def __eq__(self, other):
        return self.contigs == other.contigs

    def __hash__(self):
        return hash(str(sorted(self.contigs)))

    def __str__(self):
        return  f"{self.origin}_{self.id}  ({len(self.contigs)} contigs)"

    def overlaps_with(self, other):
        return  self.contigs & other.contigs
    
    def __and__(self, other):
        contigs = self.contigs & other.contigs
        name = f"{self.name} & {other.name}"
        origin = f"{self.origin} & {other.origin}"

        return Bin(contigs, origin, name)

    def add_length(self, length):
        self.length = length

    def intersection(self, *others):
        other_contigs = (o.contigs for o in others)
        contigs = self.contigs.intersection(*other_contigs)
        name = f"{self.name} & {' & '.join([other.name for other in others])}"
        origin = 'intersec'

        return Bin(contigs, origin, name)

    def difference(self, *others):
        other_contigs = (o.contigs for o in others)
        contigs = self.contigs.difference(*other_contigs)
        name = f"{self.name} & {' & '.join([other.name for other in others])}"
        origin = 'diff'

        return Bin(contigs, origin, name)

    def union(self, *others):
        other_contigs = (o.contigs for o in others)
        contigs = self.contigs.union(*other_contigs)
        name = f"{self.name} & {' & '.join([other.name for other in others])}"
        origin = 'union'

        return Bin(contigs, origin, name)


def get_bins_from_directory(bin_dir: str, set_name: str) -> list:

    bins = [] 
    
    for bin_fasta_file in os.listdir(bin_dir):

        bin_fasta_path = os.path.join(bin_dir, bin_fasta_file)
        bin_name = bin_fasta_file

        contigs = {name for name, _ in pyfastx.Fasta(bin_fasta_path, build_index=False)}

        bin_obj = Bin(contigs, set_name, bin_name)
        
        bins.append(bin_obj)

    return bins


def parse_bin_directories(bin_name_to_bin_dir: dict) -> dict:

    bin_name_to_bins = {}
    
    for name, bin_dir in bin_name_to_bin_dir.items():
        bin_name_to_bins[name] = get_bins_from_directory(bin_dir, name)

    return bin_name_to_bins


def get_connected_bin_graph(bin_name_to_bins):
    G = nx.Graph()

    for set1_name, set2_name in itertools.combinations(bin_name_to_bins, 2):
        set1 = bin_name_to_bins[set1_name]
        set2 = bin_name_to_bins[set2_name]
        
        # logging.debug(f"{set1_name} vs {set2_name}")
        for bin1, bin2 in itertools.product(set1, set2):
            
            if bin1.overlaps_with(bin2):
                # logging.info(f"{bin1} overlaps with {bin2}")
                G.add_edge(bin1, bin2)
    return G


def get_all_possible_combinations(clique):
    return (c for r in range(2, len(clique)+1) for c in itertools.combinations(clique, r))

def get_intersection_bins(G):
    intersect_bins = set()
    #nx.draw_shell(G, with_labels=True)
    for clique in nx.clique.find_cliques(G):
        bins_combinations = get_all_possible_combinations(clique)
        for bins in bins_combinations:
            intersec_bin = bins[0].intersection(*bins[1:])

            intersect_bins.add(intersec_bin)

    return intersect_bins



def get_difference_bins(G):
    difference_bins = set()
    #nx.draw_shell(G, with_labels=True)
    for clique in nx.clique.find_cliques(G):
        # TODO should not use combinations but another method of itertools to get all possible combination in all possible order.
        bins_combinations = get_all_possible_combinations(clique)
        for bins in bins_combinations:
            for bin_a in bins:
                bin_diff = bin_a.difference(*(b for b in bins if b != bin_a))
                difference_bins.add(bin_diff)

    return difference_bins


def get_union_bins(G):
    union_bins = set()
    #nx.draw_shell(G, with_labels=True)
    for clique in nx.clique.find_cliques(G):
        bins_combinations = get_all_possible_combinations(clique)
        for bins in bins_combinations:
            # print(f'bins {[str(b) for b in bins]}')

            for bin_a in bins:
                bin_union = bin_a.union(*(b for b in bins if b != bin_a))
                # print("DIFF", bin_union )
                
                union_bins.add(bin_union)

    return union_bins


def create_intersec_diff_bins(G):
    new_bins = set()

    #nx.draw_shell(G, with_labels=True)
    for clique in nx.clique.find_cliques(G):
        bins_combinations = get_all_possible_combinations(clique)
        for bins in bins_combinations:

            # intersection 
            intersec_bin = bins[0].intersection(*bins[1:])
            new_bins.add(intersec_bin)

            # difference 
            for bin_a in bins:
                bin_diff = bin_a.difference(*(b for b in bins if b != bin_a))
                new_bins.add(bin_diff)
        
    return new_bins

def add_bin_size(bins, contig_to_size):
    
    for bin_obj in bins:
        length = sum((contig_to_size[c] for c in bin_obj.contigs))
        bin_obj.add_length(length)

# def get_diff_bins(connected_bins_graph):

#     #nx.draw_shell(G, with_labels=True)
#     for clique in nx.clique.find_cliques(G):
#         print('=====')
#         [print(b) for b in clique]

# def get_bins_intersections(bins):
#     for bin_obj in bins:

#     return set().intersection(*sets)

def dereplicate_bin_sets(bin_sets):
    """Dereplicate bins from different bin sets to get a non redondant bin set."""

    return set().union(*bin_sets)