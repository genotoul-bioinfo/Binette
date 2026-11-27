from collections.abc import Iterable, Mapping

import pyfastx


def parse_fasta_file(fasta_file: str, index_file: str) -> pyfastx.Fasta:
    """
    Parse a FASTA file and return a pyfastx.Fasta object.

    :param fasta_file: The path to the FASTA file.

    :return: A pyfastx.Fasta object representing the parsed FASTA file.
    """
    fa = pyfastx.Fasta(fasta_file, build_index=True, index_file=index_file)
    return fa


def make_contig_index(contigs: Iterable[str]) -> dict[str, int]:
    """
    Create an index mapping for contigs.

    :param contigs: Contig names.

    :return: A dictionary mapping each contig name to its index.
    """
    contig_to_index = {contig: index for index, contig in enumerate(contigs)}
    return contig_to_index


def apply_contig_index[T](
    contig_to_index: Mapping[str, int], contig_to_info: Mapping[str, T]
) -> dict[int, T]:
    """
    Apply the contig index mapping to the contig info dictionary.

    :param contig_to_index: A dictionary mapping contig names to their corresponding index.
    :param contig_to_info: A dictionary mapping contig names to their associated information.

    :return: A dictionary mapping contig indices to their associated information.
    """
    return {contig_to_index[contig]: info for contig, info in contig_to_info.items()}
