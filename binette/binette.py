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
import pkg_resources

from binette import contig_manager, cds, diamond, bin_quality, bin_manager, io_manager as io


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


def parse_arguments():
    """Parse script arguments."""
    program_version = pkg_resources.get_distribution("Binette").version

    parser = ArgumentParser(
        description=f"Binette version={program_version}",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    # TODO add catagory to better visualize the required and the optional args
    input_arg = parser.add_mutually_exclusive_group(required=True)

    input_arg.add_argument(
        "-d",
        "--bin_dirs",
        nargs="+",
        help="list of bin folders containing each bin in a fasta file.",
    )

    input_arg.add_argument(
        "-b",
        "--contig2bin_tables",
        nargs="+",
        help="list of contig2bin table with two columns separated\
            with a tabulation: contig, bin",
    )

    parser.add_argument("-c", "--contigs", required=True, help="Contigs in fasta format.")

    parser.add_argument(
        "-m",
        "--min_completeness",
        default=10,
        type=int,
        help="Minimum completeness required for final bin selections.",
    )

    parser.add_argument("-t", "--threads", default=1, type=int, help="Number of threads.")

    parser.add_argument("-o", "--outdir", default="results", help="Output directory.")

    parser.add_argument(
        "-w",
        "--contamination_weight",
        default=5,
        type=float,
        help="Bin are scored as follow: completeness - weight * contamination. "
             "A low contamination_weight favor complete bins over low contaminated bins.",
    )

    parser.add_argument(
        "-e",
        "--extension",
        default="fasta",
        help="Extension of fasta files in bin folders "
             "(necessary when --bin_dirs is used).",
    )

    parser.add_argument(
        "--checkm2_db",
        help="Provide a path for the CheckM2 diamond database. "
        "By default the database set via <checkm2 database> is used.",
    )

    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")

    parser.add_argument("--debug", help="active debug mode", action="store_true")

    parser.add_argument("--resume", help="active resume mode", action="store_true")

    parser.add_argument("--low_mem", help="low mem mode", action="store_true")

    parser.add_argument("--version", action="version", version=program_version)

    args = parser.parse_args()
    return args


def main():
    "Orchestrate the execution of the program"

    args = parse_arguments()

    init_logging(args.verbose, args.debug)

    # Setup input parameters #

    bin_dirs = args.bin_dirs
    contig2bin_tables = args.contig2bin_tables
    contigs_fasta = args.contigs
    threads = args.threads
    outdir = args.outdir
    low_mem = args.low_mem
    contamination_weight = args.contamination_weight

    min_completeness = args.min_completeness

    # High quality threshold used just to log number of high quality bins.
    hq_max_conta = 5
    hq_min_completeness = 90

    # Temporary files #
    out_tmp_dir = os.path.join(outdir, "temporary_files")
    os.makedirs(out_tmp_dir, exist_ok=True)

    faa_file = os.path.join(out_tmp_dir, "assembly_proteins.faa")
    diamond_result_file = os.path.join(out_tmp_dir, "diamond_result.tsv")
    diamond_log = os.path.join(out_tmp_dir, "diamond_run.log")

    # Output files #
    outdir_final_bin_set = os.path.join(outdir, "final_bins")
    os.makedirs(outdir_final_bin_set, exist_ok=True)

    final_bin_report = os.path.join(outdir, "final_bins_quality_reports.tsv")

    # Flag parameters
    resume = args.resume
    debug = args.debug

    if resume and not os.path.isfile(faa_file):
        logging.error(f"Protein file {faa_file} does not exist. Resuming is not possible")
        exit(1)

    if resume and not os.path.isfile(diamond_result_file):
        logging.error(f"Diamond result file     {diamond_result_file} does not exist. Resuming is not possible")
        exit(1)

    # Loading input bin sets

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

    logging.info("Parsing contig fasta file: {contigs_fasta}")
    contigs_object = contig_manager.parse_fasta_file(contigs_fasta)
    contig_to_length = {seq.name: len(seq) for seq in contigs_object if seq.name in contigs_in_bins}

    # Predict or reuse proteins prediction and run diamond on them
    if resume:
        logging.info(f"Parsing faa file: {faa_file}.")
        contig_to_genes = cds.parse_faa_file(faa_file)
        io.check_contig_consistency(contig_to_length, contig_to_genes, contigs_fasta, faa_file)

    else:
        contigs_iterator = (s for s in contig_manager.parse_fasta_file(contigs_fasta) if s.name in contigs_in_bins)
        contig_to_genes = cds.predict(contigs_iterator, faa_file, threads)

        if not args.checkm2_db:
            # get checkm2 db stored in checkm2 install
            diamond_db_path = diamond.get_checkm2_db()
        else:
            diamond_db_path = args.checkm2_db
            
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

    # Use contig index instead of contig name to save memory
    contig_to_index, index_to_contig = contig_manager.make_contig_index(contigs_in_bins)

    contig_to_kegg_counter = contig_manager.apply_contig_index(contig_to_index, contig_to_kegg_counter)
    contig_to_genes = contig_manager.apply_contig_index(contig_to_index, contig_to_genes)
    contig_to_length = contig_manager.apply_contig_index(contig_to_index, contig_to_length)

    bin_manager.rename_bin_contigs(original_bins, contig_to_index)

    # Extract cds metadata ##

    logging.info("Compute cds metadata.")
    (
        contig_to_cds_count,
        contig_to_aa_counter,
        contig_to_aa_length,
    ) = cds.get_contig_cds_metadata(contig_to_genes, threads)

    contig_info = {
        "contig_to_cds_count": contig_to_cds_count,
        "contig_to_aa_counter": contig_to_aa_counter,
        "contig_to_aa_length": contig_to_aa_length,
        "contig_to_kegg_counter": contig_to_kegg_counter,
        "contig_to_length": contig_to_length,
    }

    logging.info("Add size and assess quality of input bins")

    # TODO paralellize
    # original_bins = bin_quality.add_bin_metrics_in_parallel(original_bins, contig_info, threads)

    bin_quality.add_bin_metrics(original_bins, contig_info, contamination_weight, threads)

    logging.info("Create intermediate bins:")
    new_bins = bin_manager.create_intermediate_bins(bin_set_name_to_bins)

    logging.info("Assess quality for supplementary intermediate bins.")
    new_bins = bin_quality.add_bin_metrics(new_bins, contig_info, contamination_weight, threads)

    logging.info("Dereplicating input bins and new bins")
    all_bins = original_bins | new_bins

    if debug:
        all_bins_for_debug = set(all_bins)
        all_bin_compo_file = os.path.join(outdir, "all_bins_quality_reports.tsv")
        
        logging.info(f"Writing all bins in {all_bin_compo_file}")
        
        # all_high_quality_bins = [b for b in all_bins_for_debug if b.contamination <= 20 and b.completeness >= 50]
        # logging.debug(f"{len(all_high_quality_bins)} bins have contamination < 20 and completeness > 50.")
        io.write_bin_info(all_bins_for_debug, all_bin_compo_file, add_contigs=True)
        

    logging.info("Select best bins")
    selected_bins = bin_manager.select_best_bins(all_bins)
    
    logging.info(f"Filtering bins: only bins with completeness >= {min_completeness} are kept")
    selected_bins = [b for b in selected_bins if b.completeness >= min_completeness]

    logging.info(f"Writing selected bins in {final_bin_report}")
    
    for b in selected_bins:
        b.contigs = {index_to_contig[c_index] for c_index in b.contigs}
    
    io.write_bin_info(selected_bins, final_bin_report)

    io.write_bins_fasta(selected_bins, contigs_fasta, outdir_final_bin_set)

    if debug:
        for sb in selected_bins:
            if sb.completeness >= hq_min_completeness and sb.contamination <= hq_max_conta:
                logging.debug(f"{sb}, {sb.completeness}, {sb.contamination}")

    hq_bins = len(
        [sb for sb in selected_bins if sb.completeness >= hq_min_completeness and sb.contamination <= hq_max_conta]
    )
    hq_bins_single = len(
        [
            sb
            for sb in selected_bins
            if sb.completeness >= hq_min_completeness and sb.contamination <= hq_max_conta and len(sb.contigs) == 1
        ]
    )
    tresholds = f"(completeness >= {hq_min_completeness} and contamination <= {hq_max_conta})"
    logging.info(f"{hq_bins}/{len(selected_bins)} selected bins have a high quality {tresholds}.")
    logging.info(
        f"{hq_bins_single}/{len(selected_bins)} selected bins have a high quality and are made of only one contig."
    )


# If this script is run from the command line then call the main function.
if __name__ == "__main__":
    # main()
    pass
