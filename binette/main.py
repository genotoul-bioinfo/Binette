#!/usr/bin/env python
"""
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Jean Mainguy, 28 nov. 2022 
License     : MIT
Maintainer  : Jean Mainguy
Portability : POSIX
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, Action, Namespace

import sys
import logging
import os

import binette
from binette import (
    contig_manager,
    cds,
    diamond,
    bin_quality,
    bin_manager,
    io_manager as io,
)
from typing import List, Dict, Optional, Set, Tuple, Union, Sequence, Any
from pathlib import Path
import pyfastx


def init_logging(verbose, debug):
    """Initialise logging."""
    if debug:
        level = logging.DEBUG
    elif verbose:
        level = logging.INFO
    else:
        level = logging.WARNING

    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s - %(message)s",
        datefmt="[%Y-%m-%d %H:%M:%S]",
    )

    logging.info("Program started")
    logging.info(
        f'command line: {" ".join(sys.argv)}',
    )


class UniqueStore(Action):
    """
    Custom argparse action to ensure an argument is provided only once.
    """

    def __call__(
        self,
        parser: ArgumentParser,
        namespace: Namespace,
        values: Union[str, Sequence[Any], None],
        option_string: Optional[str] = None,
    ) -> None:
        """
        Ensures the argument is only used once. Raises an error if the argument appears multiple times.

        :param parser: The argparse parser instance.
        :param namespace: The namespace object that will contain the parsed arguments.
        :param values: The value associated with the argument.
        :param option_string: The option string that was used to invoke this action.
        """
        # Check if the argument has already been set
        if getattr(namespace, self.dest, self.default) is not self.default:
            parser.error(
                f"Error: The argument {option_string} can only be specified once."
            )

        # Set the argument value
        setattr(namespace, self.dest, values)


def is_valid_file(parser: ArgumentParser, arg: str) -> Path:
    """
    Validates that the provided input file exists.

    :param parser: The ArgumentParser instance handling command-line arguments.
    :param arg: The path to the file provided as an argument.
    :return: A Path object representing the valid file.
    """
    path_arg = Path(arg)

    # Check if the file exists at the provided path
    if not path_arg.exists():
        parser.error(f"Error: The specified file '{arg}' does not exist.")

    return path_arg

def parse_arguments(args):
    """Parse script arguments."""

    parser = ArgumentParser(
        description=f"Binette version={binette.__version__}",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    # ------------------------
    # Input arguments
    # ------------------------
    input_group = parser.add_argument_group("Input Arguments")
    input_arg = input_group.add_mutually_exclusive_group(required=True)

    input_arg.add_argument(
        "-d",
        "--bin-dirs",
        nargs="+",
        type=lambda x: is_valid_file(parser, x),
        action=UniqueStore,
        help="List of bin folders containing each bin in a fasta file.",
    )

    input_arg.add_argument(
        "-b",
        "--contig2bin-tables",
        nargs="+",
        action=UniqueStore,
        type=lambda x: is_valid_file(parser, x),
        help="List of contig2bin tables with two columns separated "
        "by a tabulation: contig, bin.",
    )

    input_group.add_argument(
        "-c",
        "--contigs",
        required=True,
        type=lambda x: is_valid_file(parser, x),
        help="Contigs in FASTA format.",
    )

    input_group.add_argument(
        "-p",
        "--proteins",
        type=lambda x: is_valid_file(parser, x),
        help="FASTA file of predicted proteins in Prodigal format (>contigID_geneID). "
        "Skips the gene prediction step if provided.",
    )

    # ------------------------
    # Output & runtime control
    # ------------------------
    runtime_group = parser.add_argument_group("Output and Runtime Control")

    runtime_group.add_argument(
        "-o", "--outdir", default=Path("results"), type=Path, help="Output directory."
    )

    runtime_group.add_argument(
        "-t", "--threads", default=1, type=int, help="Number of threads to use."
    )

    runtime_group.add_argument(
        "--resume",
        action="store_true",
        help="Resume mode: reuse existing temporary files if possible.",
    )

    runtime_group.add_argument(
        "-v", "--verbose", help="Increase output verbosity.", action="store_true"
    )

    runtime_group.add_argument(
        "--debug", help="Activate debug mode.", action="store_true"
    )

    runtime_group.add_argument(
        "--version", action="version", version=binette.__version__
    )

    # ------------------------
    # Bin filtering & scoring
    # ------------------------
    filter_group = parser.add_argument_group("Bin Filtering and Scoring")

    filter_group.add_argument(
        "--min-completeness",
        "--min_completeness",
        default=40,
        type=int,
        help="Minimum completeness required for intermediate bin creation and final bin selection.",
    )

    filter_group.add_argument(
        "--max-contamination",
        "--max_contamination",
        default=10,
        type=int,
        help="Maximum contamination allowed for intermediate bin creation and final bin selection.",
    )

    filter_group.add_argument(
        "--min-length",
        default=200_000,
        type=int,
        help="Minimum length (bp) required for intermediate bin creation and final bin selection.",
    )

    filter_group.add_argument(
        "--max-length",
        default=10_000_000,
        type=int,
        help="Maximum length (bp) allowed for intermediate bin creation and final bin selection.",
    )

    filter_group.add_argument(
        "-w",
        "--contamination-weight",
        default=2,
        type=float,
        help="Bins are scored as: completeness - weight * contamination. "
        "A lower weight favors completeness over low contamination.",
    )

    # ------------------------
    # Advanced options
    # ------------------------
    advanced_group = parser.add_argument_group("Advanced Options")

    advanced_group.add_argument(
        "-e",
        "--fasta-extensions",
        nargs="+",
        default={".fasta", ".fa", ".fna"},
        type=str,
        help="FASTA file extensions to search for in bin directories (used with --bin-dirs).",
    )

    advanced_group.add_argument(
        "--checkm2-db",
        type=Path,
        help="Path to CheckM2 diamond database. "
        "By default the database set via <checkm2 database> is used.",
    )

    advanced_group.add_argument(
        "--low-mem", help="Enable low-memory mode for Diamond.", action="store_true"
    )

    return parser.parse_args(args)


def parse_input_files(
    bin_dirs: List[Path],
    contig2bin_tables: List[Path],
    contigs_fasta: Path,
    fasta_extensions: Set[str] = {".fasta", ".fna", ".fa"},
):
    """
    Parses input files to retrieve information related to bins and contigs.

    :param bin_dirs: List of paths to directories containing bin FASTA files.
    :param contig2bin_tables: List of paths to contig-to-bin tables.
    :param contigs_fasta: Path to the contigs FASTA file.
    :param temporary_dir: Path to the temporary directory to store intermediate files.
    :param fasta_extensions: Possible fasta extensions to look for in the bin directory.

    :return: A tuple containing:
        - List of original bins.
        - Dictionary mapping bins to lists of contigs.
        - Dictionary mapping contig names to their lengths.
    """

    if bin_dirs:
        logging.info("Parsing bin directories.")
        bin_name_to_bin_dir = io.infer_bin_set_names_from_input_paths(bin_dirs)
        bin_set_name_to_bins_info = bin_manager.parse_bin_directories(
            bin_name_to_bin_dir, fasta_extensions
        )
    else:
        logging.info("Parsing bin2contig files.")
        bin_name_to_bin_table = io.infer_bin_set_names_from_input_paths(
            contig2bin_tables
        )
        bin_set_name_to_bins_info = bin_manager.parse_contig2bin_tables(
            bin_name_to_bin_table
        )

    logging.info(f"Processing {len(bin_set_name_to_bins_info)} bin sets.")
    for bin_set_id, bins_info in bin_set_name_to_bins_info.items():
        logging.info(f" {bin_set_id} - {len(bins_info)} bins")

    contigs_in_bins = bin_manager.get_contigs_in_bin_sets(bin_set_name_to_bins_info)
    logging.info(f"Found {len(contigs_in_bins)} contigs in input bins")

    contig_to_index = contig_manager.make_contig_index(contigs_in_bins)

    contig_key_to_bin = bin_manager.make_bins_from_bins_info(
        bin_set_name_to_bins_info, contig_to_index, are_original_bins=True
    )

    # original_bins = bin_manager.dereplicate_bin_sets(bin_set_name_to_bins.values())

    logging.info(
        f"Parsing contig fasta file to retrieve lengths of contigs: {contigs_fasta}"
    )

    contigs_in_bins_set = set(contigs_in_bins)
    contig_to_length = {
        name: len(seq)
        for name, seq in pyfastx.Fastx(contigs_fasta.as_posix())
        if name in contigs_in_bins_set
    }

    logging.debug("Parsing contig fasta is done")
    # check if all contigs from input bins are present in contigs file
    unexpected_contigs = {
        contig for contig in contigs_in_bins if contig not in contig_to_length
    }

    if len(unexpected_contigs):
        raise ValueError(
            f"{len(unexpected_contigs)} contigs from the input bins were not found in the contigs file '{contigs_fasta}'. "
            f"The missing contigs are: {', '.join(unexpected_contigs)}. Please ensure all contigs from input bins are present in contig file."
        )
    logging.debug("No unexpected contigs found.")

    contig_id_to_length = {
        contig_to_index[name]: length for name, length in contig_to_length.items()
    }

    return (
        contig_key_to_bin,
        contigs_in_bins,
        contig_id_to_length,
        contig_to_index,
    )


def manage_protein_alignement(
    faa_file: Path,
    contigs_fasta: Path,
    contigs_in_bins: Set[str],
    diamond_result_file: Path,
    checkm2_db: Optional[Path],
    threads: int,
    use_existing_protein_file: bool,
    resume_diamond: bool,
    low_mem: bool,
) -> Tuple[Dict[str, int], Dict[str, List[str]]]:
    """
    Predicts or reuses proteins prediction and runs diamond on them.

    :param faa_file: The path to the .faa file.
    :param contigs_fasta: The path to the contigs FASTA file.
    :param contigs_in_bins: Dictionary mapping bin names to lists of contigs.
    :param diamond_result_file: The path to the diamond result file.
    :param checkm2_db: The path to the CheckM2 database.
    :param threads: Number of threads for parallel processing.
    :param use_existing_protein_file: Boolean indicating whether to use an existing protein file.
    :param resume_diamond: Boolean indicating whether to resume diamond alignement.
    :param low_mem: Boolean indicating whether to use low memory mode.

    :return: A tuple containing dictionaries - contig_to_kegg_counter and contig_to_genes.
    """

    # Predict or reuse proteins prediction and run diamond on them
    if use_existing_protein_file:
        logging.info(f"Parsing faa file: {faa_file}.")
        contig_to_genes = cds.parse_faa_file(faa_file.as_posix())
        io.check_contig_consistency(
            contigs_in_bins,
            contig_to_genes,
            contigs_fasta.as_posix(),
            faa_file.as_posix(),
        )

    else:
        contigs_iterator = (
            (name, seq)
            for name, seq in pyfastx.Fastx(contigs_fasta.as_posix())
            if name in contigs_in_bins
        )
        contig_to_genes = cds.predict(contigs_iterator, faa_file.as_posix(), threads)

    if not resume_diamond:
        if checkm2_db is None:
            # get checkm2 db stored in checkm2 install
            diamond_db_path = diamond.get_checkm2_db()
        elif checkm2_db.exists():
            diamond_db_path = checkm2_db.as_posix()
        else:
            raise FileNotFoundError(checkm2_db)

        diamond_log = (
            diamond_result_file.parents[0]
            / f"{diamond_result_file.stem.split('.')[0]}.log"
        )

        diamond.run(
            faa_file.as_posix(),
            diamond_result_file.as_posix(),
            diamond_db_path,
            diamond_log.as_posix(),
            threads,
            low_mem=low_mem,
        )

    logging.info("Parsing diamond results.")
    contig_to_kegg_counter = diamond.get_contig_to_kegg_id(
        diamond_result_file.as_posix()
    )

    # Check contigs from diamond vs input assembly consistency
    io.check_contig_consistency(
        contigs_in_bins,
        contig_to_kegg_counter,
        contigs_fasta.as_posix(),
        diamond_result_file.as_posix(),
    )

    return contig_to_kegg_counter, contig_to_genes


def write_bins_fasta(
    selected_bins: List[bin_manager.Bin],
    contigs_fasta: Path,
    contigs_in_bins: dict,
    outdir: Path,
):

    for b in selected_bins:
        b.contigs = {contigs_in_bins[c_index] for c_index in b.contigs}

    outdir_final_bin_set = outdir / "final_bins"

    io.write_bins_fasta(
        selected_bins, contigs_fasta, outdir_final_bin_set, contigs_in_bins
    )

    return selected_bins


def log_selected_bin_info(
    selected_bins: List[bin_manager.Bin],
    hq_min_completeness: float,
    hq_max_conta: float,
):
    """
    Log information about selected bins based on quality thresholds.

    :param selected_bins: List of Bin objects to analyze.
    :param hq_min_completeness: Minimum completeness threshold for high-quality bins.
    :param hq_max_conta: Maximum contamination threshold for high-quality bins.

    This function logs information about selected bins that meet specified quality thresholds.
    It counts the number of high-quality bins based on completeness and contamination values.
    """

    # Log completeness and contamination in debug log
    logging.debug("High quality bins:")
    for sb in selected_bins:
        if sb.is_high_quality(
            min_completeness=hq_min_completeness, max_contamination=hq_max_conta
        ):
            logging.debug(
                f"> {sb} completeness={sb.completeness}, contamination={sb.contamination}"
            )

    # Count high-quality bins and single-contig high-quality bins
    hq_bins = len(
        [
            sb
            for sb in selected_bins
            if sb.is_high_quality(
                min_completeness=hq_min_completeness, max_contamination=hq_max_conta
            )
        ]
    )

    # Log information about high-quality bins
    thresholds = (
        f"(completeness >= {hq_min_completeness} and contamination <= {hq_max_conta})"
    )
    logging.info(
        f"{hq_bins}/{len(selected_bins)} selected bins have a high quality {thresholds}."
    )


def main():
    "Orchestrate the execution of the program"

    args = parse_arguments(
        sys.argv[1:]
    )  # sys.argv is passed in order to be able to test the function parse_arguments

    init_logging(args.verbose, args.debug)

    # High quality threshold used just to log number of high quality bins.
    hq_max_conta = 5
    hq_min_completeness = 90

    write_final_fasta_bins = True

    # Temporary files #
    out_tmp_dir: Path = args.outdir / "temporary_files"
    os.makedirs(out_tmp_dir, exist_ok=True)

    use_existing_protein_file = False

    faa_file = out_tmp_dir / "assembly_proteins.faa.gz"

    diamond_result_file = out_tmp_dir / "diamond_result.tsv.gz"

    # Output files #
    final_bin_report: Path = args.outdir / "final_bins_quality_reports.tsv"
    original_bin_report_dir: Path = args.outdir / "input_bins_quality_reports"

    if args.resume:
        io.check_resume_file(faa_file, diamond_result_file)
        use_existing_protein_file = True

    (
        contig_key_to_original_bin,
        contigs_in_bins,
        contig_to_length,
        contig_to_index,
    ) = parse_input_files(
        args.bin_dirs,
        args.contig2bin_tables,
        args.contigs,
        fasta_extensions=set(args.fasta_extensions),
    )

    if args.debug:
        index_to_contig_file = args.outdir / "index_to_contig.tsv"
        logging.info(f"Writing index to contig mapping in {index_to_contig_file}")
        with open(index_to_contig_file, "w") as flout:
            flout.write("\n".join((f"{i}\t{c}" for i, c in enumerate(contigs_in_bins))))

    original_bins = list(contig_key_to_original_bin.values())

    if args.proteins and not args.resume:
        logging.info(f"Using the provided protein sequences file: {args.proteins}")
        use_existing_protein_file = True

        cds.filter_faa_file(
            contigs_in_bins,
            input_faa_file=args.proteins,
            filtered_faa_file=faa_file,
        )

    contig_name_to_kegg_counter, contig_name_to_genes = manage_protein_alignement(
        faa_file=faa_file,
        contigs_fasta=args.contigs,
        contigs_in_bins=contigs_in_bins,
        diamond_result_file=diamond_result_file,
        checkm2_db=args.checkm2_db,
        threads=args.threads,
        use_existing_protein_file=use_existing_protein_file,
        resume_diamond=args.resume,
        low_mem=args.low_mem,
    )

    contig_to_kegg_counter = contig_manager.apply_contig_index(
        contig_to_index, contig_name_to_kegg_counter
    )
    contig_to_genes = contig_manager.apply_contig_index(
        contig_to_index, contig_name_to_genes
    )

    # Extract cds metadata ##
    logging.info("Compute cds metadata.")
    contig_metadat = cds.get_contig_cds_metadata(contig_to_genes, args.threads)

    contig_metadat["contig_to_kegg_counter"] = contig_to_kegg_counter
    contig_metadat["contig_to_length"] = contig_to_length

    logging.info("Add size and assess quality of input bins")
    bin_quality.add_bin_metrics(
        original_bins, contig_metadat, args.contamination_weight, args.threads
    )
    bin_quality.add_bin_size_and_N50(original_bins, contig_to_length)

    logging.info(
        f"Writting original input bin metrics to directory: {original_bin_report_dir}"
    )
    io.write_original_bin_metrics(original_bins, original_bin_report_dir)

    logging.info("Create intermediate bins:")

    contig_lengths = bin_quality.prepare_contig_sizes(contig_to_length)

    contig_key_to_new_bin = bin_manager.create_intermediate_bins(
        contig_key_to_original_bin,
        contig_lengths=contig_lengths,
        min_comp=args.min_completeness,
        max_conta=args.max_contamination,
        min_len=args.min_length,
        max_len=args.max_length,
    )

    logging.info(f"Assess quality for {len(contig_key_to_new_bin)} intermediate bins.")

    bin_quality.add_bin_metrics(
        contig_key_to_new_bin.values(),
        contig_metadat,
        args.contamination_weight,
        args.threads,
    )

    contig_key_to_all_bin = contig_key_to_original_bin | contig_key_to_new_bin

    bin_quality.add_bin_size_and_N50(contig_key_to_all_bin.values(), contig_to_length)

    if args.debug:
        all_bin_compo_file = args.outdir / "all_bins_quality_reports.tsv"
        logging.info(f"Writing all bins in {all_bin_compo_file}")
        io.write_bin_info(
            contig_key_to_all_bin.values(), all_bin_compo_file, add_contigs=True
        )

    selected_bins = bin_manager.select_best_bins(
        contig_key_to_all_bin,
        min_completeness=args.min_completeness,
        max_contamination=args.max_contamination,
    )

    logging.info(f"Writing selected bins in {final_bin_report}")
    io.write_bin_info(selected_bins, output=final_bin_report)

    if write_final_fasta_bins:
        io.write_bins_fasta(
            selected_bins,
            args.contigs,
            outdir=args.outdir / "final_bins",
            contigs_names=contigs_in_bins,
        )

    log_selected_bin_info(selected_bins, hq_min_completeness, hq_max_conta)

    return 0
