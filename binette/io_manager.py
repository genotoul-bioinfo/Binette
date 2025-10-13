import logging
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

import pandas as pd
import pyfastx

from binette.bin_manager import Bin

logger = logging.getLogger(__name__)


def get_paths_common_prefix_suffix(
    paths: list[Path],
) -> tuple[list[str], list[str], list[str]]:
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
        common_prefix = [
            common_prefix[i]
            for i in range(common_prefix_length)
            if common_prefix[: i + 1] == part_tuple[: i + 1]
        ]
        if not common_prefix:
            break

    # Determine common suffix
    for part_tuple in parts[1:]:
        common_suffix_length = min(len(common_suffix), len(part_tuple))
        common_suffix = [
            common_suffix[-i]
            for i in range(1, common_suffix_length + 1)
            if common_suffix[-i:] == part_tuple[-i:]
        ]
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
            common_extensions = [
                common_extensions[i]
                for i in range(common_extension_length)
                if common_extensions[i] == path.suffixes[i]
            ]
            if not common_extensions:
                break

    return common_prefix, common_suffix, common_extensions


def infer_bin_set_names_from_input_paths(input_bins: list[Path]) -> dict[str, Path]:
    """
    Infer bin set names from a list of bin input directories or files.

    :param input_bins: List of input bin directories or files.
    :return: Dictionary mapping inferred bin names to their corresponding directories or files.
    """
    bin_name_to_bin_dir = {}

    common_prefix, common_suffix, common_extensions = get_paths_common_prefix_suffix(
        input_bins
    )

    for path in input_bins:
        specific_parts = path.parts[
            len(common_prefix) : len(path.parts) - len(common_suffix)
        ]

        if not common_suffix and common_extensions:
            last_specific_part = specific_parts[-1].split(".")[
                : -len(common_extensions)
            ]
            specific_parts = list(specific_parts[:-1]) + last_specific_part

        bin_set_name = "/".join(specific_parts)
        if bin_set_name == "":
            bin_set_name = path.as_posix()

        bin_name_to_bin_dir[bin_set_name] = path

    logger.debug(f"Input bins: {' '.join([path.as_posix() for path in input_bins])}")
    logger.debug(f"Common prefix to remove: {common_prefix}")
    logger.debug(f"Common suffix to remove: {common_suffix}")
    logger.debug(f"Common extension to remove: {common_suffix}")
    logger.debug(f"bin_name_to_bin_dir: {bin_name_to_bin_dir}")

    return bin_name_to_bin_dir


def write_bin_info(bins: Iterable[Bin], output: Path, add_contigs: bool = False):
    """
    Write bin information to a TSV file.

    :param bins: List of Bin objects.
    :param output: Output file path for writing the TSV.
    :param add_contigs: Flag indicating whether to include contig information.
    """

    # Define columns for the DataFrame
    columns = [
        "name",
        "origin",
        "is_original",
        "original_name",
        "completeness",
        "contamination",
        "score",
        "checkm2_model",
        "size",
        "N50",
        "coding_density",
        "contig_count",
    ]
    if add_contigs:
        columns.append("contigs")

    # Create a list of dictionaries to build the DataFrame
    data = []
    for bin_obj in sorted(
        bins, key=lambda x: (-x.score, -x.N50, -x.is_original, x.contigs_key)
    ):
        original_name = bin_obj.original_name if bin_obj.original_name else bin_obj.name
        origins = bin_obj.origin if bin_obj.is_original else {"binette"}

        bin_info = {
            "name": bin_obj.name,
            "origin": ";".join(origins),
            "is_original": bin_obj.is_original,
            "original_name": original_name,
            "completeness": bin_obj.completeness,
            "contamination": bin_obj.contamination,
            "score": round(bin_obj.score, 2),
            "checkm2_model": bin_obj.checkm2_model,
            "size": bin_obj.length,
            "N50": bin_obj.N50,
            "coding_density": round(bin_obj.coding_density, 4)
            if bin_obj.coding_density is not None
            else None,
            "contig_count": len(bin_obj.contigs),
        }

        if add_contigs:
            bin_info["contigs"] = ";".join(str(c) for c in bin_obj.contigs)

        data.append(bin_info)

    # Create pandas DataFrame and write to TSV
    df = pd.DataFrame(data, columns=columns)
    df.to_csv(output, sep="\t", index=False)


def write_bins_fasta(
    selected_bins: list[Bin],
    contigs_fasta: Path,
    outdir: Path,
    contigs_names: list[str],
    max_buffer_size: int = 50_000_000,
):
    """
    Write selected bins' contigs to separate FASTA files using pyfastx.Fastx (no index).
    Buffer entries by total character size, not just number of sequences.

    :param selected_bins: List of Bin objects with .id and .contigs.
    :param contigs_fasta: Path to the input FASTA file.
    :param outdir: Directory to save bin FASTA files.
    :param max_buffer_size: Maximum total character size to buffer before flushing.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Clear existing files for selected bins
    for sbin in selected_bins:
        out_path = outdir / f"{sbin.name}.fa"
        if out_path.exists():
            out_path.unlink()  # remove the file

    # Map contig name to bin IDs
    contig_to_bins = {}
    for sbin in selected_bins:
        for contig_id in sbin.contigs:
            contig_name = contigs_names[contig_id]
            contig_to_bins[contig_name] = sbin.name

    assert len(contig_to_bins) == sum(len(sbin.contigs) for sbin in selected_bins), (
        "Some contigs are present in multiple bins but should be unique."
    )

    buffer = defaultdict(list)
    buffer_size = 0

    def flush_buffer():
        nonlocal buffer_size
        for bin_name, seqs in buffer.items():
            if seqs:
                with open(outdir / f"{bin_name}.fa", "a") as f:
                    f.writelines(seqs)
        buffer.clear()
        buffer_size = 0

    for name, seq in pyfastx.Fastx(contigs_fasta.as_posix()):
        bin_name = contig_to_bins.get(name)
        if not bin_name:
            continue

        fasta_entry = f">{name}\n{seq}\n"
        entry_size = len(fasta_entry)

        buffer[bin_name].append(fasta_entry)

        buffer_size += entry_size
        if buffer_size >= max_buffer_size:
            flush_buffer()

    flush_buffer()


def check_contig_consistency(
    contigs_from_assembly: Iterable[str],
    contigs_from_elsewhere: Iterable[str],
    assembly_file: str,
    elsewhere_file: str,
):
    """
    Check the consistency of contig names between different sources.

    :param contigs_from_assembly: List of contig names from the assembly file.
    :param contigs_from_elsewhere: List of contig names from an external source.
    :param assembly_file: Path to the assembly file.
    :param elsewhere_file: Path to the file from an external source.
    :raises AssertionError: If inconsistencies in contig names are found.
    """
    logger.debug("Checking contig consistency")
    are_contigs_consistent = len(
        set(contigs_from_elsewhere) | set(contigs_from_assembly)
    ) <= len(set(contigs_from_assembly))

    issue_countigs = len(set(contigs_from_elsewhere) - set(contigs_from_assembly))

    message = (
        f"{issue_countigs} contigs found in file '{elsewhere_file}' "
        f"were not found in assembly_file '{assembly_file}'"
    )
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
        error_msg = (
            f"Protein file '{faa_file}' does not exist. Resuming is not possible."
        )
        logger.error(error_msg)
        raise FileNotFoundError(error_msg)

    if not diamond_result_file.exists():
        error_msg = f"Diamond result file '{diamond_result_file}' does not exist. Resuming is not possible."
        logger.error(error_msg)
        raise FileNotFoundError(error_msg)


def write_contig2bin_table(
    selected_bins: list[Bin],
    output_file: Path,
    contigs_names: list[str],
):
    """
    Write a simple TSV file mapping contig IDs to bin IDs.

    :param selected_bins: List of selected Bin objects.
    :param output_file: Path to the output TSV file.
    :param contigs_names: List of contig names where index corresponds to contig ID.
    """
    logger.info(f"Writing contig2bin table to '{output_file}'")

    # Ensure output directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        # Write contig to bin mappings
        for bin_obj in selected_bins:
            for contig_index in bin_obj.contigs:
                contig_name = contigs_names[contig_index]
                f.write(f"{contig_name}\t{bin_obj.name}\n")

    total_entries = sum(len(bin_obj.contigs) for bin_obj in selected_bins)
    logger.debug(f"Successfully wrote contig2bin table with {total_entries} entries")


def write_original_bin_metrics(original_bins: list[Bin], original_bin_report_dir: Path):
    """
    Write metrics of original input bins to a specified directory.

    This function writes the metrics for each bin set to a TSV file in the specified directory.
    Each bin set will have its own TSV file named according to its set name.

    :param original_bins: A set containing input bins
    :param original_bin_report_dir: The directory path (Path) where the bin metrics will be saved.
    """

    original_bin_report_dir.mkdir(parents=True, exist_ok=True)

    bin_set_name_to_bins = defaultdict(list)
    for bin_obj in original_bins:
        for origin in bin_obj.origin:
            bin_set_name_to_bins[origin].append(bin_obj)

    for i, (set_name, bins) in enumerate(sorted(bin_set_name_to_bins.items())):
        bins_metric_file = (
            original_bin_report_dir
            / f"input_bins_{i + 1}.{set_name.replace('/', '_')}.tsv"
        )

        logger.debug(
            f"Writing metrics for bin set '{set_name}' to file '{bins_metric_file}'"
        )
        write_bin_info(bins, bins_metric_file)

    logger.debug("Completed writing all original input bin metrics")
