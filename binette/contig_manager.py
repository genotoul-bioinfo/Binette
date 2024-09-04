import pyfastx
from typing import Dict, Iterable, Tuple, Set, Any, Union


def parse_fasta_file(fasta_file: str) -> pyfastx.Fasta:
    """
    Parse a FASTA file and return a pyfastx.Fasta object.

    :param fasta_file: The path to the FASTA file.

    :return: A pyfastx.Fasta object representing the parsed FASTA file.
    """
    fa = pyfastx.Fasta(fasta_file, build_index=True)
    return fa


def make_contig_index(contigs: Set[str]) -> Tuple[Dict[str, int], Dict[int, str]]:
    """
    Create an index mapping for contigs.

    :param contigs: A list of contig names.

    :return: A tuple containing the contig index mapping dictionaries (contig_to_index, index_to_contig).
    """
    contig_to_index = {contig: index for index, contig in enumerate(contigs)}
    index_to_contig = {index: contig for contig, index in contig_to_index.items()}
    return contig_to_index, index_to_contig


def apply_contig_index(contig_to_index: Dict[str, int], contig_to_info: Dict[str, Any]) -> Dict[int, Union[Any,Iterable[Any]]]:
    """
    Apply the contig index mapping to the contig info dictionary.

    :param contig_to_index: A dictionary mapping contig names to their corresponding index.
    :param contig_to_info: A dictionary mapping contig names to their associated information.

    :return: A dictionary mapping contig indices to their associated information.
    """
    return {contig_to_index[contig]: info for contig, info in contig_to_info.items()}
