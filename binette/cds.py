
import concurrent.futures as cf
import logging
from collections import Counter, defaultdict
from typing import Dict, List, Iterator, Tuple

import pyfastx
import pyrodigal
from tqdm import tqdm


def get_contig_from_cds_name(cds_name: str) -> str:
    """
    Extract the contig name from a CDS name.

    :param cds_name: The name of the CDS.
    :type cds_name: str
    :return: The name of the contig.
    :rtype: str
    """
    return "_".join(cds_name.split("_")[:-1])

def predict(contigs_iterator: Iterator, outfaa: str, threads: int =1) -> Dict[str, List[str]]:
    """
    Predict open reading frames with Pyrodigal.

    :param contigs_iterator: An iterator of contig sequences.
    :param outfaa: The output file path for predicted protein sequences (in FASTA format).
    :param threads: Number of CPU threads to use (default is 1).

    :return: A dictionary mapping contig names to predicted genes.
    """
    future_per_contig = {}
    orf_finder = pyrodigal.GeneFinder(meta="meta")

    logging.info(f"Predicting CDS sequences with Pyrodigal using {threads} threads.")
    with cf.ProcessPoolExecutor(max_workers=threads) as tpe:
        for seq in contigs_iterator:
            future_per_contig[seq.name] = tpe.submit(orf_finder.find_genes, seq.seq)

    contig_to_pyrodigal_genes = {contig_id: future.result() for contig_id, future in future_per_contig.items()}
    write_faa(outfaa, contig_to_pyrodigal_genes)

    contig_to_genes = {
        contig_id: [gene.translate() for gene in pyrodigal_genes]
        for contig_id, pyrodigal_genes in contig_to_pyrodigal_genes.items()
    }
    return contig_to_genes

def write_faa(outfaa: str, contig_to_genes: Dict[str, List[str]]) -> None:
    """
    Write predicted protein sequences to a FASTA file.

    :param outfaa: The output file path for predicted protein sequences (in FASTA format).
    :param contig_to_genes: A dictionary mapping contig names to predicted genes.

    """
    logging.info("Writing predicted protein sequences.")
    with open(outfaa, "w") as fl:
        for contig_id, genes in contig_to_genes.items():
            genes.write_translations(fl, contig_id)

def parse_faa_file(faa_file: str) -> Dict[str, List]:
    """
    Parse a FASTA file containing protein sequences and organize them by contig.

    :param faa_file: Path to the input FASTA file.
    :return: A dictionary mapping contig names to lists of protein sequences.
    """
    contig_to_genes = defaultdict(list)
    for name, seq in pyfastx.Fastx(faa_file):
        contig = get_contig_from_cds_name(name)
        contig_to_genes[contig].append(seq)

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

def get_contig_cds_metadata_flat(contig_to_genes: Dict[str, List[str]]) -> Tuple[Dict[str, int], Dict[str, Counter], Dict[str, int]]:
    """
    Calculate metadata for contigs, including CDS count, amino acid composition, and total amino acid length.

    :param contig_to_genes: A dictionary mapping contig names to lists of protein sequences.
    :return: A tuple containing dictionaries for CDS count, amino acid composition, and total amino acid length.
    """
    contig_to_cds_count = {contig: len(genes) for contig, genes in contig_to_genes.items()}

    contig_to_aa_counter = {contig: get_aa_composition(genes) for contig, genes in tqdm(contig_to_genes.items(), unit="contig")}
    logging.info("Calculating amino acid composition.")

    contig_to_aa_length = {contig: sum(counter.values()) for contig, counter in tqdm(contig_to_aa_counter.items(), unit="contig")}
    logging.info("Calculating total amino acid length.")

    return contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length

def get_contig_cds_metadata(contig_to_genes: Dict[str, List[str]], threads: int) -> Tuple[Dict[str, int], Dict[str, Counter], Dict[str, int]]:
    """
    Calculate metadata for contigs in parallel, including CDS count, amino acid composition, and total amino acid length.

    :param contig_to_genes: A dictionary mapping contig names to lists of protein sequences.
    :param threads: Number of CPU threads to use.
    :return: A tuple containing dictionaries for CDS count, amino acid composition, and total amino acid length.
    """
    contig_to_cds_count = {contig: len(genes) for contig, genes in contig_to_genes.items()}

    contig_to_future = {}
    logging.info(f"Collecting contig amino acid composition using {threads} threads.")
    with cf.ProcessPoolExecutor(max_workers=threads) as tpe:
        for contig, genes in tqdm(contig_to_genes.items()):
            contig_to_future[contig] = tpe.submit(get_aa_composition, genes)

    contig_to_aa_counter = {contig: future.result() for contig, future in tqdm(contig_to_future.items(), unit="contig")}
    logging.info("Calculating amino acid composition in parallel.")

    contig_to_aa_length = {contig: sum(counter.values()) for contig, counter in  tqdm(contig_to_aa_counter.items(), unit="contig")}
    logging.info("Calculating total amino acid length in parallel.")

    return contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length
