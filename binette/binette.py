#!/usr/bin/env python
"""
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Jean Mainguy, 28 nov. 2022 
License     : MIT
Maintainer  : Jean Mainguy
Portability : POSIX
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import sys
import logging
import os

import binette
from binette import contig_manager, cds, diamond, bin_quality, bin_manager, io_manager as io
from typing import List, Dict, Set, Tuple


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



def parse_arguments(args):
    """Parse script arguments."""

    parser = ArgumentParser(
        description=f"Binette version={binette.__version__}",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    # Input arguments category
    input_group = parser.add_argument_group('Input Arguments')
    input_arg = input_group.add_mutually_exclusive_group(required=True)

    input_arg.add_argument(
        "-d",
        "--bin_dirs",
        nargs="+",
        help="List of bin folders containing each bin in a fasta file.",
    )

    input_arg.add_argument(
        "-b",
        "--contig2bin_tables",
        nargs="+",
        help="List of contig2bin table with two columns separated\
            with a tabulation: contig, bin",
    )

    input_group.add_argument("-c", "--contigs", required=True, help="Contigs in fasta format.")

    # Other parameters category
    other_group = parser.add_argument_group('Other Arguments')

    other_group.add_argument(
        "-m",
        "--min_completeness",
        default=40,
        type=int,
        help="Minimum completeness required for final bin selections.",
    )

    other_group.add_argument("-t", "--threads", default=1, type=int, help="Number of threads to use.")

    other_group.add_argument("-o", "--outdir", default="results", help="Output directory.")

    other_group.add_argument(
        "-w",
        "--contamination_weight",
        default=2,
        type=float,
        help="Bin are scored as follow: completeness - weight * contamination. "
             "A low contamination_weight favor complete bins over low contaminated bins.",
    )

    other_group.add_argument(
        "--checkm2_db",
        help="Provide a path for the CheckM2 diamond database. "
        "By default the database set via <checkm2 database> is used.",
    )

    other_group.add_argument("--low_mem", help="Use low mem mode when running diamond", action="store_true")

    other_group.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")

    other_group.add_argument("--debug", help="Activate debug mode", action="store_true")

    other_group.add_argument("--resume",
                            action="store_true",
                            help="Activate resume mode. Binette will examine the 'temporary_files' directory "
                                "within the output directory and reuse any existing files if possible."
                        )


    other_group.add_argument("--version", action="version", version=binette.__version__)

    args = parser.parse_args(args)
    return args

def parse_input_files(bin_dirs: List[str], contig2bin_tables: List[str], contigs_fasta: str) -> Tuple[Dict[str, List], List, Dict[str, List], Dict[str, int]]:
    """
    Parses input files to retrieve information related to bins and contigs.

    :param bin_dirs: List of paths to directories containing bin FASTA files.
    :param contig2bin_tables: List of paths to contig-to-bin tables.
    :param contigs_fasta: Path to the contigs FASTA file.

    :return: A tuple containing:
        - Dictionary mapping bin set names to lists of bins.
        - List of original bins.
        - Dictionary mapping bins to lists of contigs.
        - Dictionary mapping contig names to their lengths.
    """

    if bin_dirs:
        logging.info("Parsing bin directories.")
        bin_name_to_bin_dir = io.infer_bin_name_from_bin_inputs(bin_dirs)
        bin_set_name_to_bins = bin_manager.parse_bin_directories(bin_name_to_bin_dir)
    else:
        logging.info("Parsing bin2contig files.")
        bin_name_to_bin_table = io.infer_bin_name_from_bin_inputs(contig2bin_tables)
        bin_set_name_to_bins = bin_manager.parse_contig2bin_tables(bin_name_to_bin_table)

    logging.info(f"{len(bin_set_name_to_bins)} bin sets processed:")
    for bin_set_id, bins in bin_set_name_to_bins.items():
        logging.info(f" {bin_set_id} - {len(bins)} bins")

    original_bins = bin_manager.dereplicate_bin_sets(bin_set_name_to_bins.values())
    contigs_in_bins = bin_manager.get_contigs_in_bins(original_bins)

    logging.info(f"Parsing contig fasta file: {contigs_fasta}")
    contigs_object = contig_manager.parse_fasta_file(contigs_fasta)
    contig_to_length = {seq.name: len(seq) for seq in contigs_object if seq.name in contigs_in_bins}

    return bin_set_name_to_bins, original_bins, contigs_in_bins, contig_to_length


def manage_protein_alignement(faa_file: str, contigs_fasta: str, contig_to_length: Dict[str, List],
                                contigs_in_bins: Dict[str, List], diamond_result_file: str,
                                checkm2_db: str, threads: int, resume: bool, low_mem: bool) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Predicts or reuses proteins prediction and runs diamond on them.
    
    :param faa_file: The path to the .faa file.
    :param contigs_fasta: The path to the contigs FASTA file.
    :param contig_to_length: Dictionary mapping contig names to their lengths.
    :param contigs_in_bins: Dictionary mapping bin names to lists of contigs.
    :param diamond_result_file: The path to the diamond result file.
    :param checkm2_db: The path to the CheckM2 database.
    :param threads: Number of threads for parallel processing.
    :param resume: Boolean indicating whether to resume the process.
    :param low_mem: Boolean indicating whether to use low memory mode.

    :return: A tuple containing dictionaries - contig_to_kegg_counter and contig_to_genes.
    """

    # Predict or reuse proteins prediction and run diamond on them
    if resume:
        logging.info(f"Parsing faa file: {faa_file}.")
        contig_to_genes = cds.parse_faa_file(faa_file)
        io.check_contig_consistency(contig_to_length, contig_to_genes, contigs_fasta, faa_file)

    else:
        contigs_iterator = (s for s in contig_manager.parse_fasta_file(contigs_fasta) if s.name in contigs_in_bins)
        contig_to_genes = cds.predict(contigs_iterator, faa_file, threads)

        if checkm2_db:
            diamond_db_path = checkm2_db
        else:
            # get checkm2 db stored in checkm2 install
            diamond_db_path = diamond.get_checkm2_db()
        
        diamond_log = f"{os.path.splitext(diamond_result_file)[0]}.log"

        diamond.run(
            faa_file,
            diamond_result_file,
            diamond_db_path,
            diamond_log,
            threads,
            low_mem=low_mem,
        )

    logging.info("Parsing diamond results.")
    contig_to_kegg_counter = diamond.get_contig_to_kegg_id(diamond_result_file)

    # Check contigs from diamond vs input assembly consistency
    io.check_contig_consistency(contig_to_length, contig_to_kegg_counter, contigs_fasta, diamond_result_file)

    return contig_to_kegg_counter, contig_to_genes


def select_bins_and_write_them(all_bins: Set[bin_manager.Bin], contigs_fasta: str, final_bin_report: str, min_completeness: float,
                               index_to_contig: dict, outdir: str, debug: bool) -> List[bin_manager.Bin]:
    """
    Selects and writes bins based on specific criteria.

    :param all_bins: Set of Bin objects.
    :param contigs_fasta: Path to the contigs FASTA file.
    :param final_bin_report: Path to write the final bin report.
    :param min_completeness: Minimum completeness threshold for bin selection.
    :param index_to_contig: Dictionary mapping indices to contig names.
    :param outdir: Output directory to save final bins and reports.
    :param debug: Debug mode flag.
    :return: Selected bins that meet the completeness threshold.
    """

    outdir_final_bin_set = os.path.join(outdir, "final_bins")
    os.makedirs(outdir_final_bin_set, exist_ok=True)

    if debug:
        all_bins_for_debug = set(all_bins)
        all_bin_compo_file = os.path.join(outdir, "all_bins_quality_reports.tsv")
        
        logging.info(f"Writing all bins in {all_bin_compo_file}")
        
        io.write_bin_info(all_bins_for_debug, all_bin_compo_file, add_contigs=True)
        
        with open(os.path.join(outdir, "index_to_contig.tsv"), 'w') as flout:
            flout.write('\n'.join((f'{i}\t{c}' for i, c in index_to_contig.items())))

    logging.info("Selecting best bins")
    selected_bins = bin_manager.select_best_bins(all_bins)

    logging.info(f"Bin Selection: {len(selected_bins)} selected bins")

    logging.info(f"Filtering bins: only bins with completeness >= {min_completeness} are kept")
    selected_bins = [b for b in selected_bins if b.completeness >= min_completeness]

    logging.info(f"Filtering bins: {len(selected_bins)} selected bins")

    logging.info(f"Writing selected bins in {final_bin_report}")
    
    for b in selected_bins:
        b.contigs = {index_to_contig[c_index] for c_index in b.contigs}
    
    io.write_bin_info(selected_bins, final_bin_report)

    io.write_bins_fasta(selected_bins, contigs_fasta, outdir_final_bin_set)

    return selected_bins



def log_selected_bin_info(selected_bins: List[bin_manager.Bin], hq_min_completeness: float, hq_max_conta: float):
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
        if sb.completeness >= hq_min_completeness and sb.contamination <= hq_max_conta:
            logging.debug(f"> {sb} completeness={sb.completeness}, contamination={sb.contamination}")

    # Count high-quality bins and single-contig high-quality bins
    hq_bins = len([sb for sb in selected_bins if sb.completeness >= hq_min_completeness and sb.contamination <= hq_max_conta])

    # Log information about high-quality bins
    thresholds = f"(completeness >= {hq_min_completeness} and contamination <= {hq_max_conta})"
    logging.info(f"{hq_bins}/{len(selected_bins)} selected bins have a high quality {thresholds}.")


def main():
    "Orchestrate the execution of the program"

    args = parse_arguments(sys.argv[1:]) # sys.argv is passed in order to be able to test the function parse_arguments

    init_logging(args.verbose, args.debug)

    # High quality threshold used just to log number of high quality bins.
    hq_max_conta = 5
    hq_min_completeness = 90

    # Temporary files #
    out_tmp_dir = os.path.join(args.outdir, "temporary_files")
    os.makedirs(out_tmp_dir, exist_ok=True)

    faa_file = os.path.join(out_tmp_dir, "assembly_proteins.faa")
    diamond_result_file = os.path.join(out_tmp_dir, "diamond_result.tsv")

    # Output files #
    final_bin_report = os.path.join(args.outdir, "final_bins_quality_reports.tsv")


    if args.resume:
        io.check_resume_file(faa_file, diamond_result_file)

    bin_set_name_to_bins, original_bins, contigs_in_bins, contig_to_length = parse_input_files(args.bin_dirs, args.contig2bin_tables, args.contigs)

    contig_to_kegg_counter, contig_to_genes = manage_protein_alignement(faa_file=faa_file, contigs_fasta=args.contigs, contig_to_length=contig_to_length,
                                contigs_in_bins=contigs_in_bins,
                                diamond_result_file=diamond_result_file, checkm2_db=args.checkm2_db,
                                threads=args.threads, resume=args.resume, low_mem=args.low_mem)
    
    # Use contig index instead of contig name to save memory
    contig_to_index, index_to_contig = contig_manager.make_contig_index(contigs_in_bins)

    contig_to_kegg_counter = contig_manager.apply_contig_index(contig_to_index, contig_to_kegg_counter)
    contig_to_genes = contig_manager.apply_contig_index(contig_to_index, contig_to_genes)
    contig_to_length = contig_manager.apply_contig_index(contig_to_index, contig_to_length)

    bin_manager.rename_bin_contigs(original_bins, contig_to_index)


    # Extract cds metadata ##
    logging.info("Compute cds metadata.")
    contig_metadat = cds.get_contig_cds_metadata(contig_to_genes, args.threads)

    contig_metadat["contig_to_kegg_counter"] = contig_to_kegg_counter
    contig_metadat["contig_to_length"] = contig_to_length
    

    logging.info("Add size and assess quality of input bins")
    bin_quality.add_bin_metrics(original_bins, contig_metadat, args.contamination_weight, args.threads)

    logging.info("Create intermediate bins:")
    new_bins = bin_manager.create_intermediate_bins(bin_set_name_to_bins)

    logging.info("Assess quality for supplementary intermediate bins.")
    new_bins = bin_quality.add_bin_metrics(new_bins, contig_metadat, args.contamination_weight, args.threads)


    logging.info("Dereplicating input bins and new bins")
    all_bins = original_bins | new_bins

    selected_bins = select_bins_and_write_them(all_bins,  args.contigs, final_bin_report, args.min_completeness, index_to_contig,  args.outdir, args.debug)

    log_selected_bin_info(selected_bins, hq_min_completeness, hq_max_conta)

    return 0