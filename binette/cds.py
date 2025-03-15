import concurrent.futures as cf
import multiprocessing.pool
import logging
from collections import Counter, defaultdict
from typing import Dict, List, Iterator, Tuple, Any, Union, Set

import pyfastx
import pyrodigal
from tqdm import tqdm
from pathlib import Path
import gzip


def get_contig_from_cds_name(cds_name: str) -> str:
    """
    Extract the contig name from a CDS name.

    :param cds_name: The name of the CDS.
    :type cds_name: str
    :return: The name of the contig.
    :rtype: str
    """
    return "_".join(cds_name.split("_")[:-1])


def predict(
    contigs_iterator: Iterator, outfaa: str, threads: int = 1
) -> Dict[str, List[str]]:
    """
    Predict open reading frames with Pyrodigal.

    :param contigs_iterator: An iterator of contig sequences.
    :param outfaa: The output file path for predicted protein sequences (in FASTA format).
    :param threads: Number of CPU threads to use (default is 1).

    :return: A dictionary mapping contig names to predicted genes.
    """
    try:
        # for version >=3 of pyrodigal
        orf_finder = pyrodigal.GeneFinder(meta="meta")  # type: ignore
    except AttributeError:
        orf_finder = pyrodigal.OrfFinder(meta="meta")  # type: ignore

    logging.info(f"Predicting cds sequences with Pyrodigal using {threads} threads.")

    with multiprocessing.pool.ThreadPool(processes=threads) as pool:
        contig_and_genes = pool.starmap(
            predict_genes,
            ((orf_finder.find_genes, name, seq) for name, seq in contigs_iterator),
        )

    write_faa(outfaa, contig_and_genes)

    contig_to_genes = {
        contig_id: [gene.translate() for gene in pyrodigal_genes]
        for contig_id, pyrodigal_genes in contig_and_genes
    }

    return contig_to_genes


def predict_genes(find_genes, name, seq) -> Tuple[str, pyrodigal.Genes]:

    return (name, find_genes(seq))


def write_faa(outfaa: str, contig_to_genes: List[Tuple[str, pyrodigal.Genes]]) -> None:
    """
    Write predicted protein sequences to a FASTA file.

    :param outfaa: The output file path for predicted protein sequences (in FASTA format).
    :param contig_to_genes: A dictionary mapping contig names to predicted genes.

    """
    logging.info("Writing predicted protein sequences.")
    with gzip.open(outfaa, "wt") as fl:
        for contig_id, genes in contig_to_genes:
            genes.write_translations(fl, contig_id)


def is_nucleic_acid(sequence: str) -> bool:
    """
    Determines whether the given sequence is a DNA or RNA sequence.

    :param sequence: The sequence to check.
    :return: True if the sequence is a DNA or RNA sequence, False otherwise.
    """
    # Define nucleotidic bases (DNA and RNA)
    nucleotidic_bases = set("ATCGNUatcgnu")

    # Check if all characters in the sequence are valid nucleotidic bases (DNA or RNA)
    if all(base in nucleotidic_bases for base in sequence):
        return True

    # If any character is invalid, return False
    return False


def parse_faa_file(faa_file: str) -> Dict[str, List[str]]:
    """
    Parse a FASTA file containing protein sequences and organize them by contig.

    :param faa_file: Path to the input FASTA file.
    :return: A dictionary mapping contig names to lists of protein sequences.
    :raises ValueError: If the file contains nucleotidic sequences instead of protein sequences.
    """
    contig_to_genes = defaultdict(list)
    checked_sequences = []

    # Iterate through the FASTA file and parse sequences
    for name, seq in pyfastx.Fastx(faa_file):
        contig = get_contig_from_cds_name(name)
        contig_to_genes[contig].append(seq)

        # Concatenate up to the first 20 sequences for validation
        if len(checked_sequences) < 20:
            checked_sequences.append(seq)

    # Concatenate all checked sequences for a more reliable nucleic acid check
    concatenated_seq = "".join(checked_sequences)

    # Check if the concatenated sequence appears to be nucleic acid
    if is_nucleic_acid(concatenated_seq):
        raise ValueError(
            f"The file '{faa_file}' appears to contain nucleotide sequences. "
            "Ensure that the file contains valid protein sequences in FASTA format."
        )

    return dict(contig_to_genes)


def get_aa_composition(genes: List[str]) -> Counter:
    """
    Calculate the amino acid composition of a list of protein sequences.

    :param genes: A list of protein sequences.
    :return: A Counter object representing the amino acid composition.
    """
    aa_counter = Counter()
    for gene in genes:
        aa_counter += Counter(gene)

    return aa_counter


def get_contig_cds_metadata_flat(
    contig_to_genes: Dict[str, List[str]]
) -> Tuple[Dict[str, int], Dict[str, Counter], Dict[str, int]]:
    """
    Calculate metadata for contigs, including CDS count, amino acid composition, and total amino acid length.

    :param contig_to_genes: A dictionary mapping contig names to lists of protein sequences.
    :return: A tuple containing dictionaries for CDS count, amino acid composition, and total amino acid length.
    """
    contig_to_cds_count = {
        contig: len(genes) for contig, genes in contig_to_genes.items()
    }

    contig_to_aa_counter = {
        contig: get_aa_composition(genes)
        for contig, genes in tqdm(contig_to_genes.items(), unit="contig")
    }
    logging.info("Calculating amino acid composition.")

    contig_to_aa_length = {
        contig: sum(counter.values())
        for contig, counter in tqdm(contig_to_aa_counter.items(), unit="contig")
    }
    logging.info("Calculating total amino acid length.")

    return contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length


def get_contig_cds_metadata(
    contig_to_genes: Dict[int, Union[Any, List[Any]]], threads: int
) -> Dict[str, Dict]:
    """
    Calculate metadata for contigs in parallel, including CDS count, amino acid composition, and total amino acid length.

    :param contig_to_genes: A dictionary mapping contig names to lists of protein sequences.
    :param threads: Number of CPU threads to use.
    :return: A tuple containing dictionaries for CDS count, amino acid composition, and total amino acid length.
    """
    contig_to_cds_count = {
        contig: len(genes) for contig, genes in contig_to_genes.items()
    }

    contig_to_future = {}
    logging.info(f"Collecting contig amino acid composition using {threads} threads.")
    with cf.ProcessPoolExecutor(max_workers=threads) as tpe:
        for contig, genes in tqdm(contig_to_genes.items()):
            contig_to_future[contig] = tpe.submit(get_aa_composition, genes)

    contig_to_aa_counter = {
        contig: future.result()
        for contig, future in tqdm(contig_to_future.items(), unit="contig")
    }
    logging.info("Calculating amino acid composition in parallel.")

    contig_to_aa_length = {
        contig: sum(counter.values())
        for contig, counter in tqdm(contig_to_aa_counter.items(), unit="contig")
    }
    logging.info("Calculating total amino acid length in parallel.")

    contig_info = {
        "contig_to_cds_count": contig_to_cds_count,
        "contig_to_aa_counter": contig_to_aa_counter,
        "contig_to_aa_length": contig_to_aa_length,
    }

    return contig_info


def filter_faa_file(
    contigs_to_keep: Set[str],
    input_faa_file: Path,
    filtered_faa_file: Path,
):
    """
    Filters a FASTA file containing protein sequences to only include sequences
    from contigs present in the provided set of contigs (`contigs_to_keep`).

    This function processes the input FASTA file, identifies protein sequences
    originating from contigs listed in `contigs_to_keep`, and writes the filtered
    sequences to a new FASTA file. The output file supports optional `.gz` compression.

    :param contigs_to_keep: A set of contig names to retain in the output FASTA file.
    :param input_faa_file: Path to the input FASTA file containing protein sequences.
    :param filtered_faa_file: Path to the output FASTA file for filtered sequences.
                              If the filename ends with `.gz`, the output will be compressed.
    """
    # Determine whether the output file should be compressed
    proper_open = gzip.open if str(filtered_faa_file).endswith(".gz") else open

    # Initialize tracking sets for metrics
    contigs_with_genes = set()
    contigs_parsed = set()

    # Process the input FASTA file and filter sequences based on contigs_to_keep
    with proper_open(filtered_faa_file, "wt") as fl:
        for name, seq in pyfastx.Fastx(input_faa_file):
            contig = get_contig_from_cds_name(name)
            contigs_parsed.add(contig)
            if contig in contigs_to_keep:
                contigs_with_genes.add(contig)
                fl.write(f">{name}\n{seq}\n")

    # Calculate metrics
    total_contigs = len(contigs_to_keep)
    contigs_with_no_genes = total_contigs - len(contigs_with_genes)
    contigs_not_in_keep_list = len(contigs_parsed - contigs_to_keep)

    # Log the computed metrics
    logging.info(f"Processing protein sequences from '{input_faa_file}'.")
    logging.info(
        f"Filtered {input_faa_file} to retain genes from {total_contigs} contigs that are included in the input bins."
    )
    logging.debug(
        f"Found {contigs_with_no_genes}/{total_contigs}  contigs ({contigs_with_no_genes / total_contigs:.2%}) with no genes."
    )
    logging.debug(
        f"{contigs_not_in_keep_list} contigs from the input FASTA file are not in the keep list."
    )
