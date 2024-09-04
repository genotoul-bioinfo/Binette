import logging
from collections import defaultdict
from pathlib import Path

import pyfastx

import itertools
import networkx as nx
from typing import List, Dict, Iterable, Tuple, Set, Mapping

class Bin:
    counter = 0

    def __init__(self, contigs: Iterable[str], origin: str, name: str, is_original:bool=False) -> None:
        """
        Initialize a Bin object.

        :param contigs: Iterable of contig names belonging to the bin.
        :param origin: Origin/source of the bin.
        :param name: Name of the bin.
        """
        Bin.counter += 1

        self.origin = {origin}
        self.name = name
        self.id = Bin.counter
        self.contigs = set(contigs)
        self.hash = hash(str(sorted(self.contigs)))

        self.length = None
        self.N50 = None

        self.completeness = None
        self.contamination = None
        self.score = None

        self.is_original = is_original

    def __eq__(self, other) -> bool:
        """
        Compare the Bin object with another object for equality.

        :param other: The object to compare with.
        :return: True if the objects are equal, False otherwise.
        """
        return self.contigs == other.contigs

    def __hash__(self) -> int:
        """
        Compute the hash value of the Bin object.

        :return: The hash value.
        """
        return self.hash

    def __str__(self) -> str:
        """
        Return a string representation of the Bin object.

        :return: The string representation of the Bin object.
        """
        return f"Bin {self.id} from {';'.join(self.origin)}  ({len(self.contigs)} contigs)"

    def overlaps_with(self, other: 'Bin') -> Set[str]:
        """
        Find the contigs that overlap between this bin and another bin.

        :param other: The other Bin object.
        :return: A set of contig names that overlap between the bins.
        """
        return self.contigs & other.contigs

    # def __and__(self, other: 'Bin') -> 'Bin':
    #     """
    #     Perform a logical AND operation between this bin and another bin.

    #     :param other: The other Bin object.
    #     :return: A new Bin object representing the intersection of the bins.
    #     """
    #     contigs = self.contigs & other.contigs
    #     name = f"{self.name} & {other.name}"
    #     origin = "intersection"

    #     return Bin(contigs, origin, name)


    def add_length(self, length: int) -> None:
        """
        Add the length attribute to the Bin object if the provided length is a positive integer.

        :param length: The length value to add.
        :return: None
        """
        if isinstance(length, int) and length > 0:
            self.length = length
        else:
            raise ValueError("Length should be a positive integer.")

    def add_N50(self, n50: int) -> None:
        """
        Add the N50 attribute to the Bin object.

        :param n50: The N50 value to add.
        :return: None
        """
        if isinstance(n50, int) and n50 >= 0:
            self.N50 = n50
        else:
            raise ValueError("N50 should be a positive integer.")
        

    def add_quality(self, completeness: float, contamination: float, contamination_weight: float) -> None:
        """
        Set the quality attributes of the bin.

        :param completeness: The completeness value.
        :param contamination: The contamination value.
        :param contamination_weight: The weight assigned to contamination in the score calculation.
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
        name = f"{self.id} & {' & '.join([str(other.id) for other in others])}"
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
        name = f"{self.id} - {' - '.join([str(other.id) for other in others])}"
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
        name = f"{self.id} | {' | '.join([str(other.id) for other in others])}"
        origin = "union"

        return Bin(contigs, origin, name)
    

    def is_complete_enough(self, min_completeness: float) -> bool:
        """
        Determine if a bin is complete enough based on completeness threshold.

        :param min_completeness: The minimum completeness required for a bin.
        
        :raises ValueError: If completeness has not been set (is None).

        :return: True if the bin meets the min_completeness threshold; False otherwise.
        """

        if self.completeness is None:
            raise ValueError(
                f"The bin '{self.name}' with ID '{self.id}' has not been evaluated for completeness or contamination, "
                "and therefore cannot be assessed."
            )

        return self.completeness >= min_completeness
    

    def is_high_quality(self, min_completeness: float, max_contamination: float) -> bool:
        """
        Determine if a bin is considered high quality based on completeness and contamination thresholds.

        :param min_completeness: The minimum completeness required for a bin to be considered high quality.
        :param max_contamination: The maximum allowed contamination for a bin to be considered high quality.
        
        :raises ValueError: If either completeness or contamination has not been set (is None).

        :return: True if the bin meets the high quality criteria; False otherwise.
        """
        if self.completeness is None or self.contamination is None:
            raise ValueError(
                f"The bin '{self.name}' with ID '{self.id}' has not been evaluated for completeness or contamination, "
                "and therefore cannot be assessed for high quality."
            )

        return self.completeness >= min_completeness and self.contamination <= max_contamination



def get_bins_from_directory(bin_dir: Path, set_name: str, fasta_extensions: Set[str]) -> List[Bin]:
    """
    Retrieves a list of Bin objects from a directory containing bin FASTA files.

    :param bin_dir: The directory path containing bin FASTA files.
    :param set_name: The name of the set the bins belong to.
    :fasta_extensions: Possible fasta extensions to look for in the bin directory.

    :return: A list of Bin objects created from the bin FASTA files.
    """
    bins = []
    fasta_extensions |= {f".{ext}" for ext in fasta_extensions if not ext.startswith(".")} # adding a dot in case given extension are lacking one
    bin_fasta_files = (fasta_file for fasta_file in bin_dir.glob("*") if set(fasta_file.suffixes) & fasta_extensions)

    for bin_fasta_path in bin_fasta_files:

        bin_name = bin_fasta_path.name

        contigs = {name for name, _ in pyfastx.Fasta(str(bin_fasta_path), build_index=False)}

        bin_obj = Bin(contigs, set_name, bin_name)

        bins.append(bin_obj)

    return bins



def parse_bin_directories(bin_name_to_bin_dir: Dict[str, Path], fasta_extensions:Set[str]) -> Dict[str, Set[Bin]]:
    """
    Parses multiple bin directories and returns a dictionary mapping bin names to a list of Bin objects.

    :param bin_name_to_bin_dir: A dictionary mapping bin names to their respective bin directories.
    :fasta_extensions: Possible fasta extensions to look for in the bin directory.

    :return: A dictionary mapping bin names to a list of Bin objects created from the bin directories.
    """
    bin_set_name_to_bins = {}

    for name, bin_dir in bin_name_to_bin_dir.items():
        bins = get_bins_from_directory(bin_dir, name, fasta_extensions)
        set_of_bins = set(bins)
        
        # Calculate the number of duplicates
        num_duplicates = len(bins) - len(set_of_bins)
        
        if num_duplicates > 0:
            logging.warning(
                f'{num_duplicates} bins with identical contig compositions detected in bin set "{name}". '
                'These bins were merged to ensure uniqueness.'
            )

        # Store the unique set of bins
        bin_set_name_to_bins[name] = set_of_bins


    return bin_set_name_to_bins

def parse_contig2bin_tables(bin_name_to_bin_tables: Dict[str, Path]) -> Dict[str, Set['Bin']]:
    """
    Parses multiple contig-to-bin tables and returns a dictionary mapping bin names to a set of unique Bin objects.

    Logs a warning if duplicate bins are detected within a bin set.

    :param bin_name_to_bin_tables: A dictionary where keys are bin set names and values are file paths or identifiers 
                                   for contig-to-bin tables. Each table is parsed to extract Bin objects.

    :return: A dictionary where keys are bin set names and values are sets of Bin objects. Duplicates are removed based 
             on contig composition.
    """
    bin_set_name_to_bins = {}

    for name, contig2bin_table in bin_name_to_bin_tables.items():
        bins = get_bins_from_contig2bin_table(contig2bin_table, name)
        set_of_bins = set(bins)
        
        # Calculate the number of duplicates
        num_duplicates = len(bins) - len(set_of_bins)
        
        if num_duplicates > 0:
            logging.warning(
                f'{num_duplicates*2} bins with identical contig compositions detected in bin set "{name}". '
                'These bins were merged to ensure uniqueness.'
            )

        # Store the unique set of bins
        bin_set_name_to_bins[name] = set_of_bins
        
    return bin_set_name_to_bins


def get_bins_from_contig2bin_table(contig2bin_table: Path, set_name: str) -> List[Bin]:
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
            contig_name = line.strip().split()[0]
            bin_name = line.strip().split("\t")[1]
            bin_name2contigs[bin_name].add(contig_name)

    bins = []
    for bin_name, contigs in bin_name2contigs.items():
        bin_obj = Bin(contigs, set_name, bin_name)
        bins.append(bin_obj)
    return bins


def from_bin_sets_to_bin_graph(bin_name_to_bin_set: Mapping[str, Iterable[Bin]]) -> nx.Graph:
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



def get_all_possible_combinations(clique: List) -> Iterable[Tuple]:
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

    :param G: A networkx Graph representing the graph of bins.
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
            bin_a = bins.pop()
            bin_union = bin_a.union(*bins)

            if bin_union.contigs:
                union_bins.add(bin_union)

    return union_bins


def select_best_bins(bins: Set[Bin]) -> List[Bin]:
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

def group_identical_bins(bins:Iterable[Bin]) -> List[List[Bin]]:
    """
    Group identical bins together

    :param bins: list of bins

    return List of list of identical bins
    """
    binhash_to_bins = defaultdict(list)

    # Collect bins by their hash values
    for bin_obj in bins:
        binhash_to_bins[bin_obj.hash].append(bin_obj)

    return list(binhash_to_bins.values())


def dereplicate_bin_sets(bin_sets: Iterable[Set['Bin']]) -> Set['Bin']:
    """
    Consolidate bins from multiple bin sets into a single set of non-redundant bins.

    Bins with the same hash are considered duplicates. For each group of duplicates,
    the origins are merged, and only one representative bin is kept.

    :param bin_sets: An iterable of sets, where each set contains `Bin` objects. These sets are merged
                     into a single set of unique bins by consolidating bins with the same hash.

    :return: A set of `Bin` objects with duplicates removed. Each `Bin` in the resulting set has
             merged origins from the bins it was consolidated with.
    """
    all_bins = (bin_obj for bins in bin_sets for bin_obj in bins)
    list_of_identical_bins = group_identical_bins(all_bins)

    dereplicated_bins = set()

    # Merge bins with the same hash
    for identical_bins in list_of_identical_bins:
        # Select the first bin as the representative
        selected_bin = identical_bins[0]
        for bin_obj in identical_bins[1:]:
            # Merge origins of all bins with the same hash
            selected_bin.origin |= bin_obj.origin

        # Add the representative bin to the result set
        dereplicated_bins.add(selected_bin)

    return dereplicated_bins

def get_contigs_in_bin_sets(bin_set_name_to_bins: Dict[str, Set[Bin]]) -> Set[str]:
    """
    Processes bin sets to check for duplicated contigs and logs detailed information about each bin set.

    :param bin_set_name_to_bins: A dictionary where keys are bin set names and values are sets of Bin objects.

    :return:  A set of contig names found in bin sets
    """
    # To track all unique contigs across bin sets
    all_contigs_in_bins = set()

    for bin_set_name, bins in bin_set_name_to_bins.items():
        list_contigs_in_bin_sets = get_contigs_in_bins(bins)

        # Count duplicates
        contig_counts = {contig: list_contigs_in_bin_sets.count(contig) for contig in list_contigs_in_bin_sets}
        duplicated_contigs = {contig: count for contig, count in contig_counts.items() if count > 1}

        if duplicated_contigs:
            logging.warning(
                f"Bin set '{bin_set_name}' contains {len(duplicated_contigs)} duplicated contigs. "
                "Details: " + ", ".join(f"{contig} (found {count} times)" for contig, count in duplicated_contigs.items())
            )

        # Unique contigs in current bin set
        unique_contigs_in_bin_set = set(list_contigs_in_bin_sets)

        # Update global contig tracker
        all_contigs_in_bins |= unique_contigs_in_bin_set

        # Log summary for the current bin set
        logging.debug(
            f"Bin set '{bin_set_name}': {len(bins)} bins, {len(unique_contigs_in_bin_set)} unique contigs."
        )

    return all_contigs_in_bins


def get_contigs_in_bins(bins: Iterable[Bin]) -> List[str]:
    """
    Retrieves all contigs present in the given list of bins.

    :param bins: A list of Bin objects.

    :return: A list of contigs present in the bins.
    """
    return [contig for b in bins for contig in b.contigs]


def rename_bin_contigs(bins: Iterable[Bin], contig_to_index: dict):
    """
    Renames the contigs in the bins based on the provided mapping.

    :param bins: A list of Bin objects.
    :param contig_to_index: A dictionary mapping old contig names to new index names.
    """
    for b in bins:
        b.contigs = {contig_to_index[contig] for contig in b.contigs}
        b.hash = hash(str(sorted(b.contigs)))

def create_intermediate_bins(bin_set_name_to_bins: Mapping[str, Iterable[Bin]]) -> Set[Bin]:
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

