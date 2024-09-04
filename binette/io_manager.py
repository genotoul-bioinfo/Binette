import logging
import pyfastx
from typing import Iterable, List, Dict, Tuple, Set
import csv

from binette.bin_manager import Bin

from pathlib import Path

def get_paths_common_prefix_suffix(paths: List[Path]) -> Tuple[List[str], List[str], List[str]]:
    """
    Determine the common prefix parts, suffix parts, and common extensions of the last part of a list of pathlib.Path objects.

    :param paths: List of pathlib.Path objects.
    :return: A tuple containing three lists:
             - The common prefix parts.
             - The common suffix parts.
             - The common extensions of the last part of the paths.
    """
    # Extract parts for all paths
    parts = [list(path.parts) for path in paths]
    
    # Find the common prefix
    if not parts:
        return [], [], []
    
    # Initialize common prefix and suffix lists
    common_prefix = list(parts[0])
    common_suffix = list(parts[0])
    # Determine common prefix
    for part_tuple in parts[1:]:
        common_prefix_length = min(len(common_prefix), len(part_tuple))
        common_prefix = [common_prefix[i] for i in range(common_prefix_length) if common_prefix[:i+1] == part_tuple[:i+1]]
        if not common_prefix:
            break

    # Determine common suffix
    for part_tuple in parts[1:]:
        common_suffix_length = min(len(common_suffix), len(part_tuple))
        common_suffix = [common_suffix[-i] for i in range(1, common_suffix_length + 1) if common_suffix[-i:] == part_tuple[-i:]]
        if not common_suffix:
            break
    if len(parts) > 1:
        common_suffix.reverse()

    # Determine common extensions of the last part of the paths
    if len(paths) == 1:
        common_extensions = paths[0].suffixes
    else:
        common_extensions = list(paths[0].suffixes)
        for path in paths[1:]:
            common_extension_length = min(len(common_extensions), len(path.suffixes))
            common_extensions = [common_extensions[i] for i in range(common_extension_length) if common_extensions[i] == path.suffixes[i]]
            if not common_extensions:
                break
    
    return common_prefix, common_suffix, common_extensions
    
def infer_bin_set_names_from_input_paths(input_bins: List[Path]) -> Dict[str, Path]:
    """
    Infer bin set names from a list of bin input directories or files.

    :param input_bins: List of input bin directories or files.
    :return: Dictionary mapping inferred bin names to their corresponding directories or files.
    """
    bin_name_to_bin_dir = {}

    common_prefix, common_suffix, common_extensions = get_paths_common_prefix_suffix(input_bins)

    for path in input_bins:

        specific_parts = path.parts[len(common_prefix):len(path.parts)-len(common_suffix)]

        if not common_suffix and common_extensions:
            last_specific_part = specific_parts[-1].split('.')[:-len(common_extensions)] 
            specific_parts = list(specific_parts[:-1]) + last_specific_part


        bin_set_name = '/'.join(specific_parts)
        if bin_set_name == "":
            bin_set_name = path.as_posix()

        bin_name_to_bin_dir[bin_set_name] = path

    logging.debug(f"Input bins: {' '.join([path.as_posix() for path in input_bins])}")
    logging.debug(f"Common prefix to remove: {common_prefix}")
    logging.debug(f"Common suffix to remove: {common_suffix}")
    logging.debug(f"Common extension to remove: {common_suffix}")
    logging.debug(f"bin_name_to_bin_dir: {bin_name_to_bin_dir}")

    return bin_name_to_bin_dir


def write_bin_info(bins: Iterable[Bin], output: Path, add_contigs: bool = False):
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
            ';'.join(bin_obj.origin),
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


def write_bins_fasta(selected_bins: List[Bin], contigs_fasta: Path, outdir: Path):
    """
    Write selected bins' contigs to separate FASTA files.

    :param selected_bins: List of Bin objects representing the selected bins.
    :param contigs_fasta: Path to the input FASTA file containing contig sequences.
    :param outdir: Output directory to save the individual bin FASTA files.
    """

    fa = pyfastx.Fasta(contigs_fasta.as_posix(), build_index=True)

    for sbin in selected_bins:
        outfile = outdir / f"bin_{sbin.id}.fa"

        with open(outfile, "w") as outfl:
            sequences = (f">{c}\n{fa[c]}" for c in sbin.contigs)
            outfl.write("\n".join(sequences) + "\n")


def check_contig_consistency(contigs_from_assembly: Iterable[str],
                             contigs_from_elsewhere: Iterable[str],
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


def check_resume_file(faa_file: Path, diamond_result_file: Path) -> None:
    """
    Check the existence of files required for resuming the process.

    :param faa_file: Path to the protein file.
    :param diamond_result_file: Path to the Diamond result file.
    :raises FileNotFoundError: If the required files don't exist for resuming.
    """

    if faa_file.exists() and diamond_result_file.exists():
        return

    if not faa_file.exists():
        error_msg = f"Protein file '{faa_file}' does not exist. Resuming is not possible."
        logging.error(error_msg)
        raise FileNotFoundError(error_msg)

    if not diamond_result_file.exists():
        error_msg = f"Diamond result file '{diamond_result_file}' does not exist. Resuming is not possible."
        logging.error(error_msg)
        raise FileNotFoundError(error_msg)


def write_original_bin_metrics(bin_set_name_to_bins: Dict[str, Set[Bin]], original_bin_report_dir: Path):
    """
    Write metrics of original input bins to a specified directory.

    This function takes a dictionary mapping bin set names to sets of bins and writes
    the metrics for each bin set to a TSV file in the specified directory. Each bin set
    will have its own TSV file named according to its set name.

    :param bin_set_name_to_bins: A dictionary where the keys are bin set names (str) and 
                                 the values are sets of Bin objects representing bins.
    :param original_bin_report_dir: The directory path (Path) where the bin metrics will be saved.
    """

    original_bin_report_dir.mkdir(parents=True, exist_ok=True)

    for i, (set_name, bins) in enumerate(sorted(bin_set_name_to_bins.items())):
        bins_metric_file = original_bin_report_dir / f"input_bins_{i + 1}.{set_name.replace('/', '_')}.tsv"
        
        logging.debug(f"Writing metrics for bin set '{set_name}' to file: {bins_metric_file}")
        write_bin_info(bins, bins_metric_file)

    logging.debug("Completed writing all original input bin metrics.")
