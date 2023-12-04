import logging
import os
import pyfastx
from typing import List, Dict
import csv

from binette.bin_manager import Bin


def infer_bin_name_from_bin_inputs(input_bins: List[str]) -> Dict[str, str]:
    """
    Infer bin names from a list of bin input directories.

    :param input_bins: List of input bin directories.
    :return: Dictionary mapping inferred bin names to their corresponding directories.
    """
    logging.debug(f"Inferring bin names from input bins:")

    commonprefix_len = len(os.path.commonprefix(input_bins))
    reversed_strings = [s[::-1] for s in input_bins]
    commonsufix_len = len(os.path.commonprefix(reversed_strings))

    bin_name_to_bin_dir = {d[commonprefix_len: len(d) - commonsufix_len]: d for d in input_bins}

    logging.debug(f"Input bins:  {' '.join(input_bins)}")
    logging.debug(f"Common prefix to remove: {os.path.commonprefix(reversed_strings)[::-1]}")
    logging.debug(f"Common suffix to remove: {os.path.commonprefix(input_bins)}")
    logging.debug(f"bin_name_to_bin_dir: {bin_name_to_bin_dir}")

    return bin_name_to_bin_dir


def write_bin_info(bins: List[Bin], output: str, add_contigs: bool = False):
    """
    Write bin information to a TSV file.

    :param bins: List of Bin objects.
    :param output: Output file path for writing the TSV.
    :param add_contigs: Flag indicating whether to include contig information.
    """

    header = ["bin_id", "origin", "name", "completeness", "contamination", "score", "size", "N50", "contig_count"]
    if add_contigs:
        header.append('contigs')

    bin_infos = []
    for bin_obj in sorted(bins, key=lambda x: (x.score, x.N50, -x.id), reverse=True):
        bin_info = [
            bin_obj.id,
            bin_obj.origin,
            bin_obj.name,
            bin_obj.completeness,
            bin_obj.contamination,
            bin_obj.score,
            bin_obj.length,
            bin_obj.N50,
            len(bin_obj.contigs),
        ]
        if add_contigs:
            bin_info.append(";".join(str(c) for c in bin_obj.contigs) if add_contigs else "")

        bin_infos.append(bin_info)

    with open(output, "w", newline="") as fl:
        writer = csv.writer(fl, delimiter="\t")
        writer.writerow(header)
        writer.writerows(bin_infos)


def write_bins_fasta(selected_bins: List[Bin], contigs_fasta: str, outdir: str):
    """
    Write selected bins' contigs to separate FASTA files.

    :param selected_bins: List of Bin objects representing the selected bins.
    :param contigs_fasta: Path to the input FASTA file containing contig sequences.
    :param outdir: Output directory to save the individual bin FASTA files.
    """

    fa = pyfastx.Fasta(contigs_fasta, build_index=True)

    for sbin in selected_bins:
        outfile = os.path.join(outdir, f"bin_{sbin.id}.fa")

        with open(outfile, "w") as outfl:
            sequences = (f">{c}\n{fa[c]}" for c in sbin.contigs)
            outfl.write("\n".join(sequences) + "\n")


def check_contig_consistency(contigs_from_assembly: List[str],
                             contigs_from_elsewhere: List[str],
                             assembly_file: str,
                             elsewhere_file: str ):
    """
    Check the consistency of contig names between different sources.

    :param contigs_from_assembly: List of contig names from the assembly file.
    :param contigs_from_elsewhere: List of contig names from an external source.
    :param assembly_file: Path to the assembly file.
    :param elsewhere_file: Path to the file from an external source.
    :raises AssertionError: If inconsistencies in contig names are found.
    """
    logging.debug("check_contig_consistency.")
    are_contigs_consistent = len(set(contigs_from_elsewhere) | set(contigs_from_assembly)) <= len(
        set(contigs_from_assembly)
    )

    issue_countigs = len(set(contigs_from_elsewhere) - set(contigs_from_assembly))
    
    message = f"{issue_countigs} contigs found in file {elsewhere_file} \
                were not found in assembly_file ({assembly_file})."
    assert are_contigs_consistent, message


def check_resume_file(faa_file: str, diamond_result_file: str) -> None:
    """
    Check the existence of files required for resuming the process.

    :param faa_file: Path to the protein file.
    :param diamond_result_file: Path to the Diamond result file.
    :raises FileNotFoundError: If the required files don't exist for resuming.
    """

    if os.path.isfile(faa_file) and os.path.isfile(diamond_result_file):
        return

    if not os.path.isfile(faa_file):
        error_msg = f"Protein file '{faa_file}' does not exist. Resuming is not possible."
        logging.error(error_msg)
        raise FileNotFoundError(error_msg)

    if not os.path.isfile(diamond_result_file):
        error_msg = f"Diamond result file '{diamond_result_file}' does not exist. Resuming is not possible."
        logging.error(error_msg)
        raise FileNotFoundError(error_msg)


