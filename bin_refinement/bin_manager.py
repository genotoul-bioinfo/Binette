"""


"""

import logging
import os
import pyfastx

import itertools
import networkx as nx


# from memory_profiler import profile

class Bin:
    counter = 0
    def __init__(self, contigs, origin, name):
        Bin.counter += 1 

        self.origin = origin
        self.name = name
        self.id = Bin.counter
        self.contigs = set(contigs)
        self.hash = hash(str(sorted(self.contigs)))

        self.length = None
        self.N50 = None
        
        self.completeness = None
        self.contamination = None
        self.score = None

    def __eq__(self, other):
        return self.contigs == other.contigs

    def __hash__(self):
        return self.hash 

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
    
    def add_N50(self, n50):
        self.N50 = n50
          

    def add_quality(self, completeness, contamination):
        self.completeness =  completeness
        self.contamination = contamination 
        self.score = completeness - 5*contamination

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
        name = f"{self.name} | {' | '.join([other.name for other in others])}"
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


def from_bin_sets_to_bin_graph(bin_name_to_bin_set):
    G = nx.Graph()

    for set1_name, set2_name in itertools.combinations(bin_name_to_bin_set, 2):
        set1 = bin_name_to_bin_set[set1_name]
        set2 = bin_name_to_bin_set[set2_name]
        
        
        for bin1, bin2 in itertools.product(set1, set2):
            
            if bin1.overlaps_with(bin2):
                G.add_edge(bin1, bin2)
    return G



def get_bin_graph(bins):
    G = nx.Graph()
    G.add_nodes_from((b.id for b in bins))     

    for i, (bin1, bin2) in enumerate(itertools.combinations(bins, 2)):

            if bin1.overlaps_with(bin2):
                # logging.info(f"{bin1} overlaps with {bin2}")
                G.add_edge(bin1.id, bin2.id, )
    return G


def get_bin_graph_with_attributes(bins, contig_to_length):
    G = nx.Graph()
    G.add_nodes_from((b.id for b in bins))     

    for i, (bin1, bin2) in enumerate(itertools.combinations(bins, 2)):
            if bin1.overlaps_with(bin2):
                
                contigs = (bin1.contigs & bin2.contigs)
                shared_length = sum((contig_to_length[c] for c in contigs))
                max_shared_length_prct = 100 - 100 * (shared_length / min((bin1.length, bin2.length)))

                # logging.info(f"{bin1} overlaps with {bin2}")
                G.add_edge(bin1.id, bin2.id, weight=max_shared_length_prct )
    return G

def get_all_possible_combinations(clique):
    return (c for r in range(2, len(clique)+1) for c in itertools.combinations(clique, r))


def get_intersection_bins(G):
    intersect_bins = set()
    #nx.draw_shell(G, with_labels=True)
    for clique in nx.clique.find_cliques(G):
        bins_combinations = get_all_possible_combinations(clique)
        for bins in bins_combinations:

            if max((b.completeness for b in bins)) < 20:
                logging.debug('completeness is not good enough to create a new bin on intersection')
                logging.debug(f"{[(str(b), b.completeness, b.contamination)  for b in bins]}")

                continue
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
                if bin_a.completeness < 20:
                    logging.debug(f'completeness of {bin_a} is not good enough to do difference... ' )
                    logging.debug(f"{[(str(b), b.completeness, b.contamination)  for b in bins]}")
                    continue
  
                if bin_diff.contigs:
                    difference_bins.add(bin_diff)

    return difference_bins


def get_union_bins(G, max_conta = 50):
    union_bins = set()
    for clique in nx.clique.find_cliques(G):
        bins_combinations = get_all_possible_combinations(clique)
        for bins in bins_combinations:
            if max((b.contamination for b in bins)) > max_conta:
                logging.debug(f"some bin are too contaminated to make a useful union bin")
                logging.debug(f"{[(str(b), b.completeness, b.contamination)  for b in bins]}")
                continue

            bins = set(bins)
            bin_a = set(bins).pop()
            bin_union = bin_a.union(*bins)
            union_bins.add(bin_union)

    return union_bins


def create_intersec_diff_bins(G):
    new_bins = set()

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
    
# @profile
# def select_best_bins(bins):
#     logging.info(f'Building bin graph from {len(bins)} ins')
#     G = get_bin_graph(bins)

#     nx.write_edgelist(G, "bin_graph_edglist")
 
#     sorted_bins = sorted(bins, key=lambda x: x.score, reverse=True)
#     selected_bins = []
#     for b in sorted_bins:
#         if b.id in G:
#             selected_bins.append(b)
#             neigbors_bins = nx.neighbors(G, b.id)
#             G.remove_nodes_from(list(neigbors_bins))
    
#     return selected_bins
    
def select_best_bins(bins):
    logging.info(f'Building no bin graph from {len(bins)} ins')

    logging.info('SORTNG')
    # sort on score, N50 and id. 
    # smaller id are prefered to select in priority original bins.
    sorted_bins = sorted(bins, key=lambda x: (x.score, x.N50, -x.id), reverse=True)

    logging.info('SELECTING')
    selected_bins = []
    for b in sorted_bins:
        if b in bins:
            overlapping_bins = {b2 for b2 in bins if b.overlaps_with(b2)}
            bins -= overlapping_bins

            selected_bins.append(b) 

    logging.info(f'SELECTING {len(selected_bins)} bins')
    return selected_bins


def dereplicate_bin_sets(bin_sets):
    """Dereplicate bins from different bin sets to get a non redondant bin set."""
    return set().union(*bin_sets)


def get_contigs_in_bins(bins):
    return set().union(*(b.contigs for b in bins))