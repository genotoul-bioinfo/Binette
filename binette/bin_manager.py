"""


"""

import logging
import os
from collections import defaultdict
import pyfastx

import itertools
import networkx as nx
from typing import List, Dict, Iterable, Tuple, Set


class Bin:
    counter = 0

    def __init__(self, contigs: Iterable[str], origin: str, name: str) -> None:
        """
        Initialize a Bin object.

        Args:
            contigs (Iterable[str]): Contig names belonging to the bin.
            origin (str): Origin/source of the bin.
            name (str): Name of the bin.

        Returns:
            None
        """
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

    def __eq__(self, other: 'Bin') -> bool:
        """
        Compare the Bin object with another object for equality.

        
        other (Any): The object to compare with.

        Returns:
            bool: True if the objects are equal, False otherwise.
        """
        return self.contigs == other.contigs

    def __hash__(self) -> int:
        """
        Compute the hash value of the Bin object.

        Returns:
            int: The hash value.
        """
        return self.hash

    def __str__(self) -> str:
        """
        Return a string representation of the Bin object.

        :return: The string representation of the Bin object.
        """
        return f"{self.origin}_{self.id}  ({len(self.contigs)} contigs)"

    def overlaps_with(self, other: 'Bin') -> Set[str]:
        """
        Find the contigs that overlap between this bin and another bin.

        :param other: The other Bin object.
        :return: A set of contig names that overlap between the bins.
        """
        return self.contigs & other.contigs

    def __and__(self, other: 'Bin') -> 'Bin':
        """
        Perform a logical AND operation between this bin and another bin.

        :param other: The other Bin object.
        :return: A new Bin object representing the intersection of the bins.
        """
        contigs = self.contigs & other.contigs
        name = f"{self.name} & {other.name}"
        origin = f"{self.origin} & {other.origin}"

        return Bin(contigs, origin, name)


    def add_length(self, length: float) -> None:
        """
        Add the length attribute to the Bin object.

        :param length: The length value to add.
        :return: None
        """
        self.length = length

    def add_N50(self, n50: float) -> None:
        """
        Add the N50 attribute to the Bin object.

        :param n50: The N50 value to add.
        :return: None
        """
        self.N50 = n50

    def add_quality(self, completeness: float, contamination: float, contamination_weight: float) -> None:
        """
        Set the quality attributes of the bin.

        :param completeness: The completeness value.
        :param contamination: The contamination value.
        :param contamination_weight: The weight assigned to contamination in the score calculation.
        :return: None
        """
        self.completeness = completeness
        self.contamination = contamination
        self.score = completeness - contamination_weight * contamination

    def intersection(self, *others: 'Bin') -> 'Bin':
        """
        Compute the intersection of the bin with other bins.

        :param others: Other bins to compute the intersection with.
        :return: A new Bin representing the intersection of the bins.
        """
        other_contigs = (o.contigs for o in others)
        contigs = self.contigs.intersection(*other_contigs)
        name = f"{self.name} & {' & '.join([other.name for other in others])}"
        origin = "intersec"

        return Bin(contigs, origin, name)

    def difference(self, *others: 'Bin') -> 'Bin':
        """
        Compute the difference between the bin and other bins.

        :param others: Other bins to compute the difference with.
        :return: A new Bin representing the difference between the bins.
        """
        other_contigs = (o.contigs for o in others)
        contigs = self.contigs.difference(*other_contigs)
        name = f"{self.name} - {' - '.join([other.name for other in others])}"
        origin = "diff"

        return Bin(contigs, origin, name)

    def union(self, *others: 'Bin') -> 'Bin':
        """
        Compute the union of the bin with other bins.

        :param others: Other bins to compute the union with.
        :return: A new Bin representing the union of the bins.
        """
        other_contigs = (o.contigs for o in others)
        contigs = self.contigs.union(*other_contigs)
        name = f"{self.name} | {' | '.join([other.name for other in others])}"
        origin = "union"

        return Bin(contigs, origin, name)




def get_bins_from_directory(bin_dir: str, set_name: str) -> List[Bin]:
    """
    Retrieves a list of Bin objects from a directory containing bin FASTA files.

    :param bin_dir: The directory path containing bin FASTA files.
    :param set_name: The name of the set the bins belong to.

    :return: A list of Bin objects created from the bin FASTA files.
    """


def parse_bin_directories(bin_name_to_bin_dir: Dict[str, str]) -> Dict[str, list]:
    """
    Parses multiple bin directories and returns a dictionary mapping bin names to a list of Bin objects.

    :param bin_name_to_bin_dir: A dictionary mapping bin names to their respective bin directories.

    :return: A dictionary mapping bin names to a list of Bin objects created from the bin directories.
    """
    bin_name_to_bins = {}

    for name, bin_dir in bin_name_to_bin_dir.items():
        bin_name_to_bins[name] = get_bins_from_directory(bin_dir, name)

    return bin_name_to_bins


def parse_contig2bin_tables(bin_name_to_bin_tables: Dict[str, str]) -> Dict[str, list]:
    """
    Parses multiple contig-to-bin tables and returns a dictionary mapping bin names to a list of Bin objects.

    :param bin_name_to_bin_tables: A dictionary mapping bin names to their respective contig-to-bin tables.

    :return: A dictionary mapping bin names to a list of Bin objects created from the contig-to-bin tables.
    """
    bin_name_to_bins = {}

    for name, contig2bin_table in bin_name_to_bin_tables.items():
        bin_name_to_bins[name] = get_bins_from_contig2bin_table(contig2bin_table, name)

    return bin_name_to_bins


def get_bins_from_contig2bin_table(contig2bin_table: str, set_name: str) -> List[Bin]:
    """
    Retrieves a list of Bin objects from a contig-to-bin table.

    :param contig2bin_table: The path to the contig-to-bin table.
    :param set_name: The name of the set the bins belong to.

    :return: A list of Bin objects created from the contig-to-bin table.
    """
    bin_name2contigs = defaultdict(set)
    with open(contig2bin_table) as fl:
        for line in fl:
            if line.startswith("#") or line.startswith("@"):
                logging.debug(f"Ignoring a line from {contig2bin_table}: {line}")
                continue
            contig_name = line.strip().split("\t")[0]
            bin_name = line.strip().split("\t")[1]
            bin_name2contigs[bin_name].add(contig_name)

    bins = []
    for bin_name, contigs in bin_name2contigs.items():
        bin_obj = Bin(contigs, set_name, bin_name)
        bins.append(bin_obj)
    return bins


def from_bin_sets_to_bin_graph(bin_name_to_bin_set: Dict[str, set]) -> nx.Graph:
    """
    Creates a bin graph from a dictionary of bin sets.

    :param bin_name_to_bin_set: A dictionary mapping bin names to their respective bin sets.

    :return: A networkx Graph representing the bin graph of overlapping bins.
    """
    G = nx.Graph()

    for set1_name, set2_name in itertools.combinations(bin_name_to_bin_set, 2):
        set1 = bin_name_to_bin_set[set1_name]
        set2 = bin_name_to_bin_set[set2_name]

        for bin1, bin2 in itertools.product(set1, set2):

            if bin1.overlaps_with(bin2):
                G.add_edge(bin1, bin2)
    return G


def get_bin_graph(bins: List[Bin]) -> nx.Graph:
    """
    Creates a bin graph from a list of Bin objects.

    :param bins: A list of Bin objects representing bins.

    :return: A networkx Graph representing the bin graph of overlapping bins.
    """
    G = nx.Graph()
    G.add_nodes_from((b.id for b in bins))

    for i, (bin1, bin2) in enumerate(itertools.combinations(bins, 2)):

        if bin1.overlaps_with(bin2):
            # logging.info(f"{bin1} overlaps with {bin2}")
            G.add_edge(
                bin1.id,
                bin2.id,
            )
    return G


def get_bin_graph_with_attributes(bins: List[Bin], contig_to_length: Dict[str, int]) -> nx.Graph:
    """
    Creates a graph from a list of Bin objects with additional attributes.

    :param bins: A list of Bin objects representing bins.
    :param contig_to_length: A dictionary mapping contig names to their lengths.

    :return: A networkx Graph representing the bin graph with attributes.
    """
    G = nx.Graph()
    G.add_nodes_from((b.id for b in bins))

    for i, (bin1, bin2) in enumerate(itertools.combinations(bins, 2)):
        if bin1.overlaps_with(bin2):

            contigs = bin1.contigs & bin2.contigs
            shared_length = sum((contig_to_length[c] for c in contigs))
            max_shared_length_prct = 100 - 100 * (shared_length / min((bin1.length, bin2.length)))

            G.add_edge(bin1.id, bin2.id, weight=max_shared_length_prct)
    return G


def get_all_possible_combinations(clique: Iterable) -> Iterable[Tuple]:
    """
    Generates all possible combinations of elements from a given clique.

    :param clique: An iterable representing a clique.

    :return: An iterable of tuples representing all possible combinations of elements from the clique.
    """
    return (c for r in range(2, len(clique) + 1) for c in itertools.combinations(clique, r))


def get_intersection_bins(G: nx.Graph) -> Set[Bin]:
    """
    Retrieves the intersection bins from a given graph.

    :param G: A networkx Graph representing the graph.

    :return: A set of Bin objects representing the intersection bins.
    """
    intersect_bins = set()

    for clique in nx.clique.find_cliques(G):
        bins_combinations = get_all_possible_combinations(clique)
        for bins in bins_combinations:
            if max((b.completeness for b in bins)) < 20:
                logging.debug("completeness is not good enough to create a new bin on intersection")
                logging.debug(f"{[(str(b), b.completeness, b.contamination)  for b in bins]}")
                continue

            intersec_bin = bins[0].intersection(*bins[1:])

            if intersec_bin.contigs:
                intersect_bins.add(intersec_bin)

    return intersect_bins


def get_difference_bins(G: nx.Graph) -> Set[Bin]:
    """
    Retrieves the difference bins from a given graph.

    :param G: A networkx Graph representing the graph.

    :return: A set of Bin objects representing the difference bins.
    """
    difference_bins = set()
    
    for clique in nx.clique.find_cliques(G):
        # TODO should not use combinations but another method of itertools
        # to get all possible combination in all possible order.
        bins_combinations = get_all_possible_combinations(clique)
        for bins in bins_combinations:

            for bin_a in bins:

                bin_diff = bin_a.difference(*(b for b in bins if b != bin_a))
                if bin_a.completeness < 20:
                    logging.debug(f"completeness of {bin_a} is not good enough to do difference... ")
                    logging.debug(f"{[(str(b), b.completeness, b.contamination)  for b in bins]}")
                    continue

                if bin_diff.contigs:
                    difference_bins.add(bin_diff)

    return difference_bins


def get_union_bins(G: nx.Graph, max_conta: int = 50) -> Set[Bin]:
    """
    Retrieves the union bins from a given graph.

    :param G: A networkx Graph representing the graph.
    :param max_conta: Maximum allowed contamination value for a bin to be included in the union.

    :return: A set of Bin objects representing the union bins.
    """
    union_bins = set()
    for clique in nx.clique.find_cliques(G):
        bins_combinations = get_all_possible_combinations(clique)
        for bins in bins_combinations:
            if max((b.contamination for b in bins)) > max_conta:
                logging.debug("Some bin are too contaminated to make a useful union bin")
                logging.debug(f"{[(str(b), b.completeness, b.contamination)  for b in bins]}")
                continue

            bins = set(bins)
            bin_a = set(bins).pop()
            bin_union = bin_a.union(*bins)

            if bin_union.contigs:
                union_bins.add(bin_union)

    return union_bins


def create_intersec_diff_bins(G: nx.Graph) -> Set[Bin]:
    """
    Creates intersection and difference bins from a given graph.

    :param G: A networkx Graph representing the graph.

    :return: A set of Bin objects representing the intersection and difference bins.
    """
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

def select_best_bins(bins: List[Bin]) -> List[Bin]:
    """
    Selects the best bins from a list of bins based on their scores, N50 values, and IDs.

    :param bins: A list of Bin objects.

    :return: A list of selected Bin objects.
    """

    logging.info("Sorting bins")
    # Sort on score, N50, and ID. Smaller ID values are preferred to select original bins first.
    sorted_bins = sorted(bins, key=lambda x: (x.score, x.N50, -x.id), reverse=True)

    logging.info("Selecting bins")
    selected_bins = []
    for b in sorted_bins:
        if b in bins:
            overlapping_bins = {b2 for b2 in bins if b.overlaps_with(b2)}
            bins -= overlapping_bins

            selected_bins.append(b)

    logging.info(f"Selected {len(selected_bins)} bins")
    return selected_bins


def dereplicate_bin_sets(bin_sets):
    """
    Dereplicates bins from different bin sets to obtain a non-redundant bin set.

    :param bin_sets: A list of bin sets.

    :return: A set of non-redundant bins.
    """
    return set().union(*bin_sets)


def get_contigs_in_bins(bins: List[Bin]) -> Set[str]:
    """
    Retrieves all contigs present in the given list of bins.

    :param bins: A list of Bin objects.

    :return: A set of contigs present in the bins.
    """
    return set().union(*(b.contigs for b in bins))


def rename_bin_contigs(bins: List[Bin], contig_to_index: dict):
    """
    Renames the contigs in the bins based on the provided mapping.

    :param bins: A list of Bin objects.
    :param contig_to_index: A dictionary mapping old contig names to new index names.
    """
    for b in bins:
        b.contigs = {contig_to_index[contig] for contig in b.contigs}


def create_intermediate_bins(bin_set_name_to_bins: Dict[str, Set[Bin]]) -> Set[Bin]:
    """
    Creates intermediate bins from a dictionary of bin sets.

    :param bin_set_name_to_bins: A dictionary mapping bin set names to corresponding bins.

    :return: A set of intermediate bins created from intersections, differences, and unions.
    """
    logging.info("Making bin graph...")
    connected_bins_graph = from_bin_sets_to_bin_graph(bin_set_name_to_bins)

    logging.info("Creating intersection bins...")
    intersection_bins = get_intersection_bins(connected_bins_graph)
    logging.info(f"{len(intersection_bins)} bins created on intersections.")

    logging.info("Creating difference bins...")
    difference_bins = get_difference_bins(connected_bins_graph)
    logging.info(f"{len(difference_bins)} bins created based on symmetric difference.")

    logging.info("Creating union bins...")
    union_bins = get_union_bins(connected_bins_graph)
    logging.info(f"{len(union_bins)} bins created on unions.")

    return difference_bins | intersection_bins | union_bins

