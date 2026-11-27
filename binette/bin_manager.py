import hashlib
import itertools
import logging
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping
from functools import cached_property
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pyfastx
from pyroaring import BitMap

logger = logging.getLogger(__name__)


class Bin:
    CHECKM2_MODELS = (
        "Neural Network (Specific Model)",
        "Gradient Boost (General Model)",
    )

    def __init__(
        self,
        contigs: BitMap,
        origin: str | None = None,
        name: str | None = None,
        is_original: bool = False,
    ) -> None:
        """
        Initialize a Bin object.

        :param contigs: Iterable of contig names belonging to the bin.
        :param origin: Origin/source of the bin.
        :param name: Name of the bin.
        """

        if not isinstance(contigs, BitMap):
            raise TypeError("Contigs should be a BitMap object.")

        self.origin: set[str] = set() if origin is None else {origin}

        self.name = name

        self.is_original = is_original

        self.contigs = contigs

        self.length: int | None = None
        self.N50: int | None = None

        self.completeness: float | None = None
        self.contamination: float | None = None
        self.score: float | None = None
        self.coding_density: float | None = None
        self.original_name: str | None = None
        self._checkm2_model_index: int | None = None

    @cached_property
    def contigs_key(self) -> bytes:
        """
        Serialize the contigs for easier comparison.
        """
        return self.contigs.serialize()

    def __eq__(self, other: object) -> bool:
        """
        Compare the Bin object with another object for equality.

        :param other: The object to compare with.
        :return: True if the objects are equal, False otherwise.
        """
        if not isinstance(other, Bin):
            return NotImplemented
        return self.contigs_key == other.contigs_key

    def __str__(self) -> str:
        """
        Return a string representation of the Bin object.

        :return: The string representation of the Bin object.
        """
        return f"Bin {self.name} from {';'.join(self.origin)}  ({len(self.contigs)} contigs)"

    def overlaps_with(self, other: "Bin") -> BitMap:
        """
        Find the contigs that overlap between this bin and another bin.

        :param other: The other Bin object.
        :return: The contigs shared by the two bins.
        """
        return self.contigs & other.contigs

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

    def add_quality(
        self, completeness: float, contamination: float, contamination_weight: float
    ) -> None:
        """
        Set the quality attributes of the bin.

        :param completeness: The completeness value.
        :param contamination: The contamination value.
        :param contamination_weight: The weight assigned to contamination in the score calculation.
        """
        self.completeness = completeness
        self.contamination = contamination
        self.score = completeness - contamination_weight * contamination

    def add_model(self, model: str) -> None:
        """
        Add a CheckM2 model to the bin.

        :param model: The model name to add.
        :raises ValueError: If the model name is not recognized.
        """

        try:
            self._checkm2_model_index = self.CHECKM2_MODELS.index(model)
        except ValueError as exc:
            raise ValueError(
                f"Unknown model '{model}' attempted to be added to bin '{self.name}'. "
                f"Valid models are: {', '.join(self.CHECKM2_MODELS)}"
            ) from exc

    @cached_property
    def checkm2_model(self) -> str | None:
        """
        Get the CheckM2 model for the bin.
        """
        if self._checkm2_model_index is not None:
            return self.CHECKM2_MODELS[self._checkm2_model_index]
        return None

    def contig_intersection(self, *others: "Bin") -> BitMap:
        """
        Compute the intersection of the bin with other bins.

        :param others: Other bins to compute the intersection with.
        """
        return self.contigs.intersection(*(o.contigs for o in others))

    def contig_difference(self, *others: "Bin") -> BitMap:
        """
        Compute the difference between the bin and other bins.

        :param others: Other bins to compute the difference with.
        """
        return self.contigs.difference(*(o.contigs for o in others))

    def contig_union(self, *others: "Bin") -> BitMap:
        """
        Compute the union of the bin with other bins.

        :param others: Other bins to compute the union with.
        """
        return self.contigs.union(*(o.contigs for o in others))

    def is_high_quality(
        self, min_completeness: float, max_contamination: float
    ) -> bool:
        """
        Determine if a bin is considered high quality based on completeness and contamination thresholds.

        :param min_completeness: The minimum completeness required for a bin to be considered high quality.
        :param max_contamination: The maximum allowed contamination for a bin to be considered high quality.

        :raises ValueError: If either completeness or contamination has not been set (is None).

        :return: True if the bin meets the high quality criteria; False otherwise.
        """
        if self.completeness is None or self.contamination is None:
            raise ValueError(
                f"The bin '{self.name}' with ID '{self.name}' has not been evaluated for completeness or contamination, "
                "and therefore cannot be assessed for high quality."
            )

        return (
            self.completeness >= min_completeness
            and self.contamination <= max_contamination
        )

    def add_coding_density(
        self, contig_to_coding_length: Mapping[int, int | None]
    ) -> float | None:
        """
        Calculate the coding density of the bin.

        :param contig_to_coding_length: A dictionary mapping contig IDs to their total coding lengths.

        :return: The coding density of the bin, or None if the length is not set or is zero.
        """
        if self.length is None or self.length == 0:
            return None

        total_coding = sum(
            contig_to_coding_length.get(contig_id) or 0 for contig_id in self.contigs
        )

        self.coding_density = total_coding / self.length


def make_bins_from_bins_info(
    bin_set_name_to_bins_info: dict[str, list[dict[str, Any]]],
    contig_to_index: dict[str, int],
    are_original_bins: bool,
) -> dict[bytes, Bin]:
    """
    Create Bin objects from the provided bin information.

    :param bin_set_name_to_bins_info: A dictionary mapping bin set names to their bin information.
    :param contig_to_index: A mapping of contig names to their indices.
    :param are_original_bins: A boolean indicating whether the bins are original.

    :return: A dictionary mapping serialized contig bitmaps to their corresponding Bin objects.
    """

    contig_key_to_bin: dict[bytes, Bin] = {}

    for set_name, bins_info in bin_set_name_to_bins_info.items():
        for bin_info in bins_info:
            bitmap_contigs = BitMap(
                contig_to_index[contig] for contig in bin_info["contigs"]
            )

            bin_obj = Bin(
                contigs=bitmap_contigs,
                origin=set_name,
                name=bin_info["bin_name"],
                is_original=are_original_bins,
            )
            if bin_obj.contigs_key not in contig_key_to_bin:
                contig_key_to_bin[bin_obj.contigs_key] = bin_obj
            else:
                # Same contigs as a bin already seen: keep that bin and record
                # that this bin set also produced it.
                contig_key_to_bin[bin_obj.contigs_key].origin.add(set_name)

    return contig_key_to_bin


def get_bins_from_directory(
    bin_dir: Path, set_name: str, fasta_extensions: set[str]
) -> list[Bin]:
    """
    Retrieves a list of Bin objects from a directory containing bin FASTA files.

    :param bin_dir: The directory path containing bin FASTA files.
    :param set_name: The name of the set the bins belong to.
    :fasta_extensions: Possible fasta extensions to look for in the bin directory.

    :return: A list of Bin objects created from the bin FASTA files.
    """
    bins = []
    fasta_extensions |= {
        f".{ext}" for ext in fasta_extensions if not ext.startswith(".")
    }  # adding a dot in case given extension are lacking one
    bin_fasta_files = (
        fasta_file
        for fasta_file in bin_dir.glob("*")
        if set(fasta_file.suffixes) & fasta_extensions
    )

    for bin_fasta_path in bin_fasta_files:
        bin_name = bin_fasta_path.with_suffix("").name

        contigs = {name for name, _ in pyfastx.Fastx(str(bin_fasta_path))}

        bin_dict = {"contigs": contigs, "set_name": set_name, "bin_name": bin_name}

        bins.append(bin_dict)

    return bins


def parse_bin_directories(
    bin_name_to_bin_dir: dict[str, Path], fasta_extensions: set[str]
) -> dict[str, list[dict[str, Any]]]:
    """
    Parses multiple bin directories and returns a dictionary mapping bin set names to a list of bin info dicts.

    Bins with identical contigs are merged later, in make_bins_from_bins_info.

    :param bin_name_to_bin_dir: A dictionary mapping bin names to their respective bin directories.
    :fasta_extensions: Possible fasta extensions to look for in the bin directory.

    :return: A dictionary mapping bin names to a list of dict created from the bin directories.
    """
    bin_set_name_to_bins_info = {}

    for name, bin_dir in sorted(bin_name_to_bin_dir.items()):
        bins = get_bins_from_directory(bin_dir, name, fasta_extensions)
        bin_set_name_to_bins_info[name] = bins

    return bin_set_name_to_bins_info


def parse_contig2bin_tables(
    bin_name_to_bin_tables: dict[str, Path],
) -> dict[str, list[dict[str, Any]]]:
    """
    Parses multiple contig-to-bin tables and returns a dictionary mapping bin set names to a list of bin info dicts.

    Bins with identical contigs are merged later, in make_bins_from_bins_info.

    :param bin_name_to_bin_tables: A dictionary where keys are bin set names and values are file paths or identifiers
                                   for contig-to-bin tables.

    :return: A dictionary mapping bin set names to a list of dict created from the contig-to-bin tables.
    """
    bin_set_name_to_bins_info = {}

    for name, contig2bin_table in sorted(bin_name_to_bin_tables.items()):
        bins = get_bins_from_contig2bin_table(contig2bin_table, name)
        bin_set_name_to_bins_info[name] = bins

    return bin_set_name_to_bins_info


def get_bins_from_contig2bin_table(
    contig2bin_table: Path, set_name: str
) -> list[dict[str, Any]]:
    """
    Retrieves a list of Bin objects from a contig-to-bin table.

    :param contig2bin_table: The path to the contig-to-bin table.
    :param set_name: The name of the set the bins belong to.

    :return: A list of Bin info in dict created from the contig-to-bin table.
    """
    bin_name2contigs = defaultdict(set)
    with open(contig2bin_table) as fl:
        for line in fl:
            if line.startswith(("#", "@")):
                logger.debug(f"Ignoring a line from {contig2bin_table}: {line}")
                continue
            contig_name = line.strip().split()[0]
            bin_name = line.strip().split("\t")[1]
            bin_name2contigs[bin_name].add(contig_name)

    bins = []
    for bin_name, contigs in bin_name2contigs.items():
        bin_dict = {"contigs": contigs, "set_name": set_name, "bin_name": bin_name}
        bins.append(bin_dict)
    return bins


def from_bins_to_bin_graph(bins: Iterable[Bin]) -> nx.Graph:
    """
    Creates a bin graph made of overlapping gram a set of bins.

    :param bins: a set of bins

    :return: A networkx Graph representing the bin graph of overlapping bins.
    """
    G = nx.Graph()

    for bin1, bin2 in itertools.combinations(bins, 2):
        if bin1.overlaps_with(bin2):
            G.add_edge(bin1.contigs_key, bin2.contigs_key)
    return G


def iter_cliques(graph: nx.Graph) -> Iterator[tuple]:
    """
    Yield every clique of at least two nodes of a graph, exactly once.

    These are all the subsets of at least two nodes of the maximal cliques. Going
    through the maximal cliques would yield a subset once for every maximal clique
    containing it.

    :param graph: Graph with sortable nodes.

    :return: Iterator of cliques, as tuples of nodes in sorted order.
    """
    # Each clique is built from its smallest node, adding nodes in increasing order.
    larger_neighbors = {
        node: {neighbor for neighbor in graph[node] if neighbor > node}
        for node in graph
    }

    def extend(clique: tuple, candidates: set) -> Iterator[tuple]:
        # candidates: nodes larger than the last node of the clique and connected to
        # all of its nodes.
        for node in sorted(candidates):
            larger_clique = (*clique, node)
            yield larger_clique
            yield from extend(larger_clique, candidates & larger_neighbors[node])

    for node in sorted(graph):
        yield from extend((node,), larger_neighbors[node])


def sort_bins_by_quality(bins: Iterable[Bin]) -> list[Bin]:
    """
    Sort bins from best to worst: highest score, then highest N50, then input bins
    before intermediate bins, then by serialized contig set.

    The serialized contig sets are only computed for bins that tie on the other
    criteria, and are not cached on the bins: sorting millions of bins doesn't keep a
    serialized copy of each contig set.

    :param bins: Bins to sort.

    :return: Sorted list of bins.
    """

    def quality(b: Bin) -> tuple[float, int, int]:
        assert b.score is not None, f"Bin {b.name} has no score"
        assert b.N50 is not None, f"Bin {b.name} has no N50"
        return (-b.score, -b.N50, -b.is_original)

    sorted_bins = []
    for _, tied_bins in itertools.groupby(sorted(bins, key=quality), key=quality):
        tied_bins = list(tied_bins)
        if len(tied_bins) > 1:
            tied_bins.sort(key=lambda b: b.contigs.serialize())
        sorted_bins.extend(tied_bins)
    return sorted_bins


def select_best_bins(
    all_bins: Iterable[Bin],
    min_completeness: float,
    max_contamination: float,
    prefix: str = "binette",
) -> list[Bin]:
    """
    Select the best non-overlapping bins based on score, N50, and ID.

    :param all_bins: Candidate bins, input and intermediate.
    :param min_completeness: Minimum completeness threshold for a bin to be considered.
    :param max_contamination: Maximum contamination threshold for a bin to be considered.
    :param prefix: Prefix to use for naming selected bins.

    """
    logger.info("Selecting best bins")

    logger.info(
        f"Filtering bins: only bins with completeness >= {min_completeness} "
        f"and contamination <= {max_contamination}"
    )
    good_enough_bins = (
        b for b in all_bins if b.is_high_quality(min_completeness, max_contamination)
    )

    logger.info("Sorting bins")
    sorted_bins = sort_bins_by_quality(good_enough_bins)

    logger.info("Selecting bins")
    # A bin is selected when it shares no contig with the bins selected before it.
    # Tracking the contigs already used is enough: no per-contig index of bins is
    # needed, which would take a lot of memory with many candidate bins.
    selected_bins = []
    used_contigs = BitMap()

    for bin_obj in sorted_bins:
        if bin_obj.contigs.isdisjoint(used_contigs):
            selected_bins.append(bin_obj)
            used_contigs |= bin_obj.contigs

    logger.info(f"Selected {len(selected_bins)} bins")

    for i, selected_bin in enumerate(selected_bins, start=1):
        if not selected_bin.origin:
            selected_bin.origin = {"binette"}
        if selected_bin.name is not None:
            selected_bin.original_name = selected_bin.name
        selected_bin.name = f"{prefix}_bin{i}"
    return selected_bins


def get_contigs_in_bin_sets(
    bin_set_name_to_bins: dict[str, list[dict[str, Any]]],
) -> list[str]:
    """
    Processes bin sets to check for duplicated contigs and logs detailed information about each bin set.

    :param bin_set_name_to_bins: A dictionary where keys are bin set names and values are lists of bin info
        dictionaries (with a "contigs" key).

    :return: The contig names found in bin sets.
    """
    # To track all unique contigs across bin sets
    all_contigs_in_bins = set()

    for bin_set_name, bins_info in bin_set_name_to_bins.items():
        list_contigs_in_bin_sets = [
            contig for bin_info in bins_info for contig in bin_info["contigs"]
        ]

        contig_counts = Counter(list_contigs_in_bin_sets)
        duplicated_contigs = {
            contig: count for contig, count in contig_counts.items() if count > 1
        }

        if duplicated_contigs:
            logger.warning(
                f"Bin set '{bin_set_name}' contains {len(duplicated_contigs)} duplicated contigs. "
                "Details: "
                + ", ".join(
                    f"{contig} (found {count} times)"
                    for contig, count in duplicated_contigs.items()
                )
            )

        # Unique contigs in current bin set
        unique_contigs_in_bin_set = set(list_contigs_in_bin_sets)

        # Update global contig tracker
        all_contigs_in_bins |= unique_contigs_in_bin_set

        # Log summary for the current bin set
        logger.debug(
            f"Bin set '{bin_set_name}': {len(bins_info)} bins, {len(unique_contigs_in_bin_set)} unique contigs."
        )

    return list(all_contigs_in_bins)


def sum_contig_lengths(bm_contigs: BitMap, contig_lengths: np.ndarray) -> int:
    """
    Sum the lengths of a set of contigs.

    :param bm_contigs: Contig indices.
    :param contig_lengths: Contig lengths indexed by contig index.

    :return: Total length of the contigs.
    """
    return int(contig_lengths[np.asarray(bm_contigs.to_array())].sum())


def iter_intermediate_contig_sets(
    contig_key_to_initial_bin: dict[bytes, Bin],
    contig_lengths: np.ndarray,
    min_comp: float,
    max_conta: float,
    min_len: int,
    max_len: int,
) -> Iterator[BitMap]:
    """
    Yield the contig sets of intermediate bins made from overlapping input bins.

    For every group of input bins that all overlap each other (a clique of the bin
    graph), the intersection, the differences and the union of their contigs are new
    candidate bins. Intersections and differences only use bins with completeness >=
    ``min_comp`` and length >= ``min_len``; unions only use bins with contamination <=
    ``max_conta`` and length <= ``max_len``. Candidates shorter than ``min_len`` or
    longer than ``max_len`` are discarded.

    Each contig set is yielded once, and never one of an input bin. Sets are yielded
    as they are created, so they don't all need to be held in memory.

    :param contig_key_to_initial_bin: Input bins by contig key.
    :param contig_lengths: Contig lengths indexed by contig index.
    :param min_comp: Minimum completeness of bins used for intersections and differences.
    :param max_conta: Maximum contamination of bins used for unions.
    :param min_len: Minimum length of bins used and created.
    :param max_len: Maximum length of bins used and created.

    :return: Iterator of contig sets.
    """
    logger.info("Making bin graph")
    connected_bins_graph = from_bins_to_bin_graph(contig_key_to_initial_bin.values())

    logger.info("Creating union, difference, and intersection bins")
    logger.debug(f"{min_comp} min completeness for intersection and difference bins")
    logger.debug(f"{max_conta} max contamination for intersection and difference bins")
    logger.debug(f"{min_len} min length for intersection and difference bins")
    logger.debug(f"{max_len} max length for intersection and difference bins")
    logger.info(
        f"Intermediate bins filtered by minimum length of {min_len} and maximum length of {max_len}."
    )

    created_count: Counter[str] = Counter()
    discarded_count: Counter[str] = Counter()
    # Contig sets already created or discarded. They are only looked up: a 16-byte
    # digest of their key is kept instead of the key itself.
    seen_digests: set[bytes] = set()

    def is_new_bin(contigs: BitMap, operation: str) -> bool:
        if not contigs:
            return False
        contig_key = contigs.serialize()
        if contig_key in contig_key_to_initial_bin:
            return False
        digest = hashlib.blake2b(contig_key, digest_size=16).digest()
        if digest in seen_digests:
            return False
        seen_digests.add(digest)
        if min_len <= sum_contig_lengths(contigs, contig_lengths) <= max_len:
            created_count[operation] += 1
            return True
        discarded_count[operation] += 1
        return False

    def usable_for_intersection_and_difference(b: Bin) -> bool:
        assert b.completeness is not None
        assert b.length is not None
        return b.completeness >= min_comp and b.length >= min_len

    def usable_for_union(b: Bin) -> bool:
        assert b.contamination is not None
        assert b.length is not None
        return b.contamination <= max_conta and b.length <= max_len

    for bin_contig_keys in iter_cliques(connected_bins_graph):
        bins = [contig_key_to_initial_bin[ck] for ck in bin_contig_keys]

        if all(usable_for_intersection_and_difference(b) for b in bins):
            contigs = bins[0].contig_intersection(*bins[1:])
            if is_new_bin(contigs, "intersection"):
                yield contigs

        for bin_a in bins:
            if usable_for_intersection_and_difference(bin_a):
                contigs = bin_a.contig_difference(*(b for b in bins if b != bin_a))
                if is_new_bin(contigs, "difference"):
                    yield contigs

        if all(usable_for_union(b) for b in bins):
            contigs = bins[0].contig_union(*bins[1:])
            if is_new_bin(contigs, "union"):
                yield contigs

    for operation, label in (
        ("intersection", "Intersection"),
        ("difference", "Symmetric Difference"),
        ("union", "Union"),
    ):
        logger.info(
            f"{label}: {created_count[operation]} bins created, "
            f"{discarded_count[operation]} discarded due to size constraints."
        )
    logger.info(
        f"{created_count.total()} new bins created from "
        f"{len(contig_key_to_initial_bin)} input bins."
    )


def create_intermediate_bins(
    contig_key_to_initial_bin: dict[bytes, Bin],
    contig_lengths: np.ndarray,
    min_comp: float,
    max_conta: float,
    min_len: int,
    max_len: int,
) -> dict[bytes, Bin]:
    """
    Create all intermediate bins at once (see iter_intermediate_contig_sets).

    :return: New bins by contig key.
    """
    return {
        contigs.serialize(): Bin(contigs, is_original=False)
        for contigs in iter_intermediate_contig_sets(
            contig_key_to_initial_bin,
            contig_lengths,
            min_comp=min_comp,
            max_conta=max_conta,
            min_len=min_len,
            max_len=max_len,
        )
    }
