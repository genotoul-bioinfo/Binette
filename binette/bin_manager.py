import itertools
import logging
from collections import Counter, defaultdict
from collections.abc import Iterable
from functools import cached_property
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pyfastx
from pyroaring import BitMap
from rich.progress import Progress

logger = logging.getLogger(__name__)


class Bin:
    CHECKM2_MODELS = (
        "Neural Network (Specific Model)",
        "Gradient Boost (General Model)",
    )

    def __init__(
        self,
        contigs: BitMap,
        origin: set[str] | None = None,
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

        if origin is None:
            self.origin = set()
        else:
            self.origin = {origin}

        self.name = name

        self.is_original = is_original

        self.contigs = contigs

        self.length = None
        self.N50 = None

        self.completeness = None
        self.contamination = None
        self.score = None
        self.coding_density = None
        self.original_name = None
        self._checkm2_model_index = None

    @cached_property
    def contigs_key(self):
        """
        Serialize the contigs for easier comparison.
        """
        return self.contigs.serialize()

    def __eq__(self, other: "Bin") -> bool:
        """
        Compare the Bin object with another object for equality.

        :param other: The object to compare with.
        :return: True if the objects are equal, False otherwise.
        """
        return self.contigs_key == other.contigs_key

    def __str__(self) -> str:
        """
        Return a string representation of the Bin object.

        :return: The string representation of the Bin object.
        """
        return f"Bin {self.name} from {';'.join(self.origin)}  ({len(self.contigs)} contigs)"

    def overlaps_with(self, other: "Bin") -> set[str]:
        """
        Find the contigs that overlap between this bin and another bin.

        :param other: The other Bin object.
        :return: A set of contig names that overlap between the bins.
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
    def checkm2_model(self):
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
        self, contig_to_coding_length: dict[int, int]
    ) -> float | None:
        """
        Calculate the coding density of the bin.

        :param contig_to_coding_length: A dictionary mapping contig IDs to their total coding lengths.

        :return: The coding density of the bin, or None if the length is not set or is zero.
        """
        if self.length is None or self.length == 0:
            return None

        total_coding = sum(
            contig_to_coding_length.get(contig_id, 0) for contig_id in self.contigs
        )

        self.coding_density = total_coding / self.length


def make_bins_from_bins_info(
    bin_set_name_to_bins_info: dict[str, list[dict[str, Any]]],
    contig_to_index: dict[str, int],
    are_original_bins: bool,
):
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
                bin_obj.origin.add(set_name)

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
    Parses multiple bin directories and returns a dictionary mapping bin names to a list of Bin objects.

    :param bin_name_to_bin_dir: A dictionary mapping bin names to their respective bin directories.
    :fasta_extensions: Possible fasta extensions to look for in the bin directory.

    :return: A dictionary mapping bin names to a list of dict created from the bin directories.
    """
    bin_set_name_to_bins_info = {}

    for name, bin_dir in sorted(bin_name_to_bin_dir.items()):
        bins = get_bins_from_directory(bin_dir, name, fasta_extensions)

        # TODO: redo this check
        # set_of_bins = set(bins)

        # # Calculate the number of duplicates
        # num_duplicates = len(bins) - len(set_of_bins)

        # if num_duplicates > 0:
        #     logger.warning(
        #         f'{num_duplicates} bins with identical contig compositions detected in bin set "{name}". '
        #         "These bins were merged to ensure uniqueness."
        #     )

        # Store the unique set of bins
        bin_set_name_to_bins_info[name] = bins

    return bin_set_name_to_bins_info


def parse_contig2bin_tables(
    bin_name_to_bin_tables: dict[str, Path],
) -> dict[str, list[dict[str, Any]]]:
    """
    Parses multiple contig-to-bin tables and returns a dictionary mapping bin names to a set of unique Bin objects.

    Logs a warning if duplicate bins are detected within a bin set.

    :param bin_name_to_bin_tables: A dictionary where keys are bin set names and values are file paths or identifiers
                                   for contig-to-bin tables. Each table is parsed to extract Bin objects.

    :return: A dictionary where keys are bin set names and values are sets of Bin objects. Duplicates are removed based
             on contig composition.
    """
    bin_set_name_to_bins_info = {}

    for name, contig2bin_table in sorted(bin_name_to_bin_tables.items()):
        bins = get_bins_from_contig2bin_table(contig2bin_table, name)

        # TODO: redo this check
        # set_of_bins = set(bins)

        # # Calculate the number of duplicates
        # num_duplicates = len(bins) - len(set_of_bins)

        # if num_duplicates > 0:
        #     logger.warning(
        #         f'{num_duplicates*2} bins with identical contig compositions detected in bin set "{name}". '
        #         "These bins were merged to ensure uniqueness."
        #     )

        # Store the unique set of bins
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
            if line.startswith("#") or line.startswith("@"):
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


def get_all_possible_combinations(clique: list) -> Iterable[tuple]:
    """
    Generates all possible combinations of elements from a given clique.

    :param clique: An iterable representing a clique.

    :return: An iterable of tuples representing all possible combinations of elements from the clique.
    """
    return (
        c for r in range(2, len(clique) + 1) for c in itertools.combinations(clique, r)
    )


def build_contig_index(bins_dict: dict[bytes, Bin]) -> dict[int, set[bytes]]:
    """
    Build an inverted index: contig_id -> set of contigs_key of bins containing it.
    :param bins_dict: Mapping from contigs_key -> Bin.
    :return: Inverted index (contig_id -> set of contigs_key).
    """
    contig_to_bins = defaultdict(set)
    for key, bin_obj in bins_dict.items():
        for contig in bin_obj.contigs:
            contig_to_bins[contig].add(key)
    return contig_to_bins


def remove_bins_from_index(
    bin_keys: set[bytes],
    bins_dict: dict[bytes, Bin],
    contig_to_bins: dict[int, set[bytes]],
) -> None:
    """
    Remove a set of bins from the inverted index.

    :param bin_keys: The contig_keys of bins to remove.
    :param bins_dict: Mapping from contig_key -> Bin.
    :param contig_to_bins: Inverted index (contig_id -> set of contig_keys).
    """
    for key in bin_keys:
        ob = bins_dict.get(key)
        if ob is not None:
            for c in ob.contigs:
                contig_to_bins[c].discard(key)


def select_best_bins(
    bins_dict: dict[bytes, Bin],
    min_completeness: float,
    max_contamination: float,
    prefix: str = "binette",
) -> list[Bin]:
    """
    Select the best non-overlapping bins based on score, N50, and ID.

    :param bins_dict: Mapping from contig_key -> Bin.
    :param min_completeness: Minimum completeness threshold for a bin to be considered.
    :param max_contamination: Maximum contamination threshold for a bin to be considered.
    :param prefix: Prefix to use for naming selected bins.

    """
    logger.info("Selecting best bins")

    logger.info(
        f"Filtering bins: only bins with completeness >= {min_completeness} "
        f"and contamination <= {max_contamination}"
    )
    good_enough_bins = {
        k: b
        for k, b in bins_dict.items()
        if b.completeness >= min_completeness and b.contamination <= max_contamination
    }

    logger.info("Sorting bins")
    sorted_bin_keys = sorted(
        good_enough_bins,
        key=lambda k: (
            -good_enough_bins[k].score,
            -good_enough_bins[k].N50,
            -good_enough_bins[k].is_original,
            k,  # contigs_key itself is sortable (bytes)
        ),
    )

    logger.info("Building contig index")
    contig_to_bin_keys = build_contig_index(good_enough_bins)

    logger.info("Selecting bins")
    selected_bins = []
    discarded_keys = set()

    for bin_key in sorted_bin_keys:
        if bin_key in discarded_keys:
            continue

        bin_obj = good_enough_bins[bin_key]
        selected_bins.append(bin_obj)

        # Gather overlapping bins via inverted index
        overlapping_bin_keys = set()
        for contig in bin_obj.contigs:
            overlapping_bin_keys |= contig_to_bin_keys[contig]

        # Discard them
        discarded_keys |= overlapping_bin_keys

        # Remove discarded bins from index to shrink future lookups
        remove_bins_from_index(
            overlapping_bin_keys, good_enough_bins, contig_to_bin_keys
        )

    logger.info(f"Selected {len(selected_bins)} bins")

    for i, selected_bin in enumerate(selected_bins, start=1):
        if not selected_bin.origin:
            selected_bin.origin = {"binette"}
        if selected_bin.name is not None:
            selected_bin.original_name = selected_bin.name
        selected_bin.name = f"{prefix}_bin{i}"
    return selected_bins


def get_contigs_in_bin_sets(bin_set_name_to_bins: dict[str, set[Bin]]) -> list[str]:
    """
    Processes bin sets to check for duplicated contigs and logs detailed information about each bin set.

    :param bin_set_name_to_bins: A dictionary where keys are bin set names and values are sets of Bin objects.

    :return:  A set of contig names found in bin sets
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


def sum_contig_lengths(
    bm_contigs: BitMap,
    contig_lengths: np.ndarray,
    cache: dict[bytes, int] | None = None,
    key: bytes | None = None,
):
    if cache is None:
        cache = {}
    if key is None:
        key = bm_contigs.serialize()
    if key not in cache:
        cache[key] = int(contig_lengths[np.fromiter(bm_contigs, dtype=np.int32)].sum())
    return cache[key]


def create_intermediate_bins(
    contig_key_to_initial_bin: dict[bytes, Bin],
    contig_lengths: np.ndarray,
    min_comp: float,
    max_conta: float,
    min_len: int,
    max_len: int,
    disable_progress_bar: bool = False,
) -> dict[bytes, Bin]:
    """
    Creates intermediate bins from a dictionary of bin sets.

    :param original_bins: Set of input bins.

    :return: A set of intermediate bins created from intersections, differences, and unions.
    """
    bin_length_cache = {}

    logger.info("Making bin graph")
    connected_bins_graph = from_bins_to_bin_graph(contig_key_to_initial_bin.values())

    cliques_of_bins = sorted(
        [sorted(clique) for clique in nx.clique.find_cliques(connected_bins_graph)]
    )

    logger.info("Creating union, difference, and intersection bins")
    logger.debug(f"{min_comp} min completeness for intersection and difference bins")
    logger.debug(f"{max_conta} max contamination for intersection and difference bins")
    logger.debug(f"{min_len} min length for intersection and difference bins")
    logger.debug(f"{max_len} max length for intersection and difference bins")
    logger.info(
        f"Intermediate bins filtered by minimum length of {min_len} and maximum length of {max_len}."
    )
    intersec_count = 0
    union_count = 0
    diff_count = 0

    intersec_size_discarded_count = 0
    diff_size_discarded_count = 0
    union_size_discarded_count = 0

    contig_key_to_new_contigs_set = {}
    discarded_contig_set_keys = set()
    with Progress(disable=disable_progress_bar) as progress:
        task = progress.add_task(
            f"Processing {len(cliques_of_bins)} cliques of bins",
            total=len(cliques_of_bins),
        )
        for clique in cliques_of_bins:
            progress.update(task, advance=1)
            bins_combinations = get_all_possible_combinations(clique)

            for bin_contig_keys in bins_combinations:
                bins = [contig_key_to_initial_bin[ck] for ck in bin_contig_keys]

                if all(
                    b.completeness >= min_comp and b.length >= min_len for b in bins
                ):
                    intersec_contigs = bins[0].contig_intersection(*bins[1:])

                    if intersec_contigs:
                        contig_key = intersec_contigs.serialize()

                        if (
                            contig_key not in contig_key_to_initial_bin
                            and contig_key not in contig_key_to_new_contigs_set
                            and contig_key not in discarded_contig_set_keys
                        ):
                            contigs_length = sum_contig_lengths(
                                intersec_contigs,
                                contig_lengths,
                                cache=bin_length_cache,
                                key=contig_key,
                            )

                            if contigs_length >= min_len and contigs_length <= max_len:
                                contig_key_to_new_contigs_set[contig_key] = (
                                    intersec_contigs
                                )
                                intersec_count += 1
                            else:
                                discarded_contig_set_keys.add(contig_key)
                                intersec_size_discarded_count += 1

                for bin_a in bins:
                    if bin_a.completeness >= min_comp and bin_a.length >= min_len:
                        diff_contigs = bin_a.contig_difference(
                            *(b for b in bins if b != bin_a)
                        )

                        if diff_contigs:
                            contig_key = diff_contigs.serialize()

                            if (
                                contig_key not in contig_key_to_initial_bin
                                and contig_key not in contig_key_to_new_contigs_set
                                and contig_key not in discarded_contig_set_keys
                            ):
                                contigs_length = sum_contig_lengths(
                                    diff_contigs,
                                    contig_lengths,
                                    cache=bin_length_cache,
                                    key=contig_key,
                                )

                                if (
                                    contigs_length >= min_len
                                    and contigs_length <= max_len
                                ):
                                    contig_key_to_new_contigs_set[contig_key] = (
                                        diff_contigs
                                    )
                                    diff_count += 1
                                else:
                                    discarded_contig_set_keys.add(contig_key)
                                    diff_size_discarded_count += 1

                if all(
                    b.contamination <= max_conta and b.length <= max_len for b in bins
                ):
                    union_contigs = bins[0].contig_union(*bins[1:])
                    if union_contigs:
                        contig_key = union_contigs.serialize()
                        if (
                            contig_key not in contig_key_to_initial_bin
                            and contig_key not in contig_key_to_new_contigs_set
                            and contig_key not in discarded_contig_set_keys
                        ):
                            contigs_length = sum_contig_lengths(
                                union_contigs,
                                contig_lengths,
                                cache=bin_length_cache,
                                key=contig_key,
                            )
                            if contigs_length >= min_len and contigs_length <= max_len:
                                contig_key_to_new_contigs_set[contig_key] = (
                                    union_contigs
                                )
                                union_count += 1
                            else:
                                discarded_contig_set_keys.add(contig_key)
                                union_size_discarded_count += 1

    logger.info(
        f"Intersection: {intersec_count} bins created, {intersec_size_discarded_count} discarded due to size constraints."
    )

    logger.info(
        f"Symmetric Difference: {diff_count} bins created, {diff_size_discarded_count} discarded due to size constraints."
    )

    logger.info(
        f"Union: {union_count} bins created, {union_size_discarded_count} discarded due to size constraints."
    )

    contig_key_to_new_bin: dict[bytes, Bin] = {
        contig_key: Bin(contigs, is_original=False)
        for contig_key, contigs in contig_key_to_new_contigs_set.items()
    }

    logger.info(
        f"{len(contig_key_to_new_bin)} new bins created from {len(contig_key_to_initial_bin)} input bins."
    )

    return contig_key_to_new_bin
