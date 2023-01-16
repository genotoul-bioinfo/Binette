#!/usr/bin/env python
'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Jean Mainguy, 28 nov. 2022 
License     : MIT 
Maintainer  : jean.mainguy@inrae.fr 
Portability : POSIX


'''

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import sys
import logging
import os

import contig_manager
import cds
import diamond
import bin_quality
import pyfastx
import bin_manager

import pkg_resources

PROGRAM_NAME = "Binette"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def init_logging(verbose, debug):
    '''Initialise logging.'''
    if debug:
        level = logging.DEBUG
    elif verbose:
        level = logging.INFO
    else:
        level = logging.WARNING


    logging.basicConfig(level=level,
                        format='%(asctime)s %(levelname)s - %(message)s',
                        datefmt="[%Y-%m-%d %H:%M:%S]")

    logging.info('Program started')
    logging.info(f'command line: {" ".join(sys.argv)}', )



def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(description="...",
                            formatter_class=ArgumentDefaultsHelpFormatter)

    input_arg = parser.add_mutually_exclusive_group(required=True)

    input_arg.add_argument("-d", "--bin_dirs", nargs='+', 
                        help="list of bin folders containing each bin in a fasta file.")

    input_arg.add_argument("-b", "--contig2bin_tables", nargs='+',
                        help="list of contig2bin table with two columns separated with a tabulation: contig, bin")                   
    
    parser.add_argument("-c", "--contigs", required=True, 
                        help="Contigs in fasta format.")
    parser.add_argument("-m", "--min_completeness", 
                        default=10, type=int, help="Minimum completeness required for final bin selections.")
    parser.add_argument("-t", "--threads", 
                        default=1, type=int, help="Number of threads.")

    parser.add_argument("-o", "--outdir", default='results',
                        help="Output directory.")
    parser.add_argument("-w", "--contamination_weigth", default=5, type=float,
                        help="Bin are scored as follow: completeness - weigth * contamination. A low contamination_weigth favor complete bins over low contaminated bins.")
                        

    parser.add_argument("-e", "--extension", default='fasta', 
                        help="Extension of fasta files in bin folders (necessary when --bin_dirs is used).")

    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")
    parser.add_argument("--debug", help="active debug mode",
                    action="store_true")
    parser.add_argument("--resume", help="active resume mode",
                        action="store_true")
    parser.add_argument("--low_mem", help="low mem mode",
                        action="store_true")
    parser.add_argument('--version',
                        action='version',
                        version=PROGRAM_VERSION)
                        
    args = parser.parse_args()
    return args

def infer_bin_name_from_bin_inputs(input_bins):
    # remove common prefix of bin dir to get nicer label
    commonprefix_len = len(os.path.commonprefix(input_bins))
    reversed_strings = [''.join(reversed(s)) for s in input_bins]
    commonsufix_len = len(os.path.commonprefix(reversed_strings))

    bin_name_to_bin_dir = {d[commonprefix_len:len(d)-commonsufix_len]:d for d in input_bins}

    logging.debug(f"input bins : {input_bins}")
    logging.debug(f"commonprefix {os.path.commonprefix(input_bins)}")
    logging.debug(f"commonsuffix {os.path.commonprefix(reversed_strings)}")
    logging.debug(f"bin_name_to_bin_dir {bin_name_to_bin_dir}")

    return bin_name_to_bin_dir

def write_bin_info(bins, output):

    header = ['bin_id', 'origin', "name", 'completeness', 'contamination', "score", 'size', 'N50', 'contig_count']
    with open(output, "w") as fl:
        fl.write('\t'.join(header)+"\n")
        for bin_obj in bins:

            line = [bin_obj.id, bin_obj.origin, bin_obj.name, 
                    bin_obj.completeness, bin_obj.contamination, bin_obj.score, 
                    bin_obj.length, bin_obj.N50,
                     len(bin_obj.contigs), ]

            fl.write("\t".join((str(e) for e in line)) + '\n')

def write_bin_info_debug(bins, output):

    header = ['origin', "name", "id", 'completeness', 'contamination', "score", 'size', 'N50', 'contig_count', "contigs"]

    with open(output, "w") as fl:
        fl.write('\t'.join(header)+"\n")
        for bin_obj in bins:

            line = [bin_obj.origin, bin_obj.name, bin_obj.id,
                    bin_obj.completeness, bin_obj.contamination, bin_obj.score, 
                    bin_obj.length, bin_obj.N50,
                     len(bin_obj.contigs), 
                     ";".join((str(c) for c in bin_obj.contigs))]

            fl.write("\t".join((str(e) for e in line)) + '\n')

def write_bins_fasta(selected_bins, contigs_fasta, outdir):
    fa = pyfastx.Fasta(contigs_fasta, build_index=True)

    for sbin in selected_bins:
        outfile = os.path.join(outdir, f"bin_{sbin.id}.fa")

        with open(outfile, 'w') as outfl:
            sequences = (f">{c}\n{fa[c]}" for c in sbin.contigs)
            outfl.write('\n'.join(sequences)+"\n")
            

def check_contig_consistency(contigs_from_assembly, contigs_from_elsewhere, assembly_file, elsewhere_file):
    logging.debug('check_contig_consistency.')
    are_contigs_consistent = len(set(contigs_from_elsewhere) | set(contigs_from_assembly)) <= len(set(contigs_from_assembly))

    message = f"{len(set(contigs_from_elsewhere) - set(contigs_from_assembly))} contigs found in file {elsewhere_file} were not found in assembly_file ({assembly_file})."
    assert are_contigs_consistent, message

def main():
    "Orchestrate the execution of the program"

    args = parse_arguments()

    init_logging(args.verbose, args.debug)

    ### Setup input parameters ###

    bin_dirs = args.bin_dirs
    contig2bin_tables = args.contig2bin_tables
    contigs_fasta = args.contigs
    threads = args.threads
    outdir = args.outdir
    low_mem = args.low_mem
    contamination_weigth = args.contamination_weigth

    min_completeness = args.min_completeness
    
    # High quality threshold used just to log number of high quality bins. 
    hq_max_conta = 5
    hq_min_completeness = 90

    ## Temporary files ##
    out_tmp_dir = os.path.join(outdir, 'temporary_files')
    os.makedirs(out_tmp_dir, exist_ok=True)

    faa_file = os.path.join(out_tmp_dir, 'assembly_proteins.faa')
    diamond_result_file = os.path.join(out_tmp_dir, 'diamond_result.tsv' ) 
    diamond_log  = os.path.join(out_tmp_dir, 'diamond_run.log' ) 

    ## Output files ##
    outdir_final_bin_set = os.path.join(outdir, 'final_bins')
    os.makedirs(outdir_final_bin_set, exist_ok=True)

    final_bin_report = os.path.join(outdir, 'final_bins_quality_reports.tsv')

    ## Flag parameters ##
    resume = args.resume
    debug = args.debug

    if resume and not os.path.isfile(faa_file):
        logging.error(f'Protein file {faa_file} does not exist. Resuming is not possible')
        exit(1)

    if resume and not os.path.isfile(diamond_result_file):
        logging.error(f'Diamond result file {diamond_result_file} does not exist. Resuming is not possible')
        exit(1) 

    ### Loading input bin sets ####

    if bin_dirs: 
        logging.info('Parsing bin directories.')
        bin_name_to_bin_dir = infer_bin_name_from_bin_inputs(bin_dirs)
        bin_set_name_to_bins = bin_manager.parse_bin_directories(bin_name_to_bin_dir)
    else:
        logging.info('Parsing bin2contig files.')
        bin_name_to_bin_table = infer_bin_name_from_bin_inputs(contig2bin_tables)
        bin_set_name_to_bins = bin_manager.parse_contig2bin_tables(bin_name_to_bin_table)

    logging.info(f'{len(bin_set_name_to_bins)} bin sets processed:')
    for bin_set_id, bins in bin_set_name_to_bins.items():
        logging.info(f' {bin_set_id} - {len(bins)} ')

    

    original_bins = bin_manager.dereplicate_bin_sets(bin_set_name_to_bins.values()) 
    contigs_in_bins = bin_manager.get_contigs_in_bins(original_bins) 

    logging.info('Parse fasta file of assembly.')
    contigs_object = contig_manager.parse_fasta_file(contigs_fasta)
    contig_to_length = {seq.name:len(seq) for seq in contigs_object if seq.name in contigs_in_bins}

    ### Predict or reuse proteins prediction and run diamond on them ####
    if resume:
        logging.info(f'Parsing faa file: {faa_file}.')
        contig_to_genes = cds.parse_faa_file(faa_file)
        check_contig_consistency(contig_to_length, contig_to_genes, contigs_fasta, faa_file)

    else:
        contigs_iterator = (s for s in contig_manager.parse_fasta_file(contigs_fasta) if s.name in contigs_in_bins)
        contig_to_genes = cds.predict(contigs_iterator, faa_file, threads)

        diamond_db_path = diamond.get_checkm2_db()
        diamond.run(faa_file, diamond_result_file, diamond_db_path, diamond_log, threads, low_mem=low_mem)

    logging.info('Parsing diamond results.')
    contig_to_kegg_counter = diamond.get_contig_to_kegg_id(diamond_result_file)
    ## Check contigs from diamond vs input assembly consistency
    check_contig_consistency(contig_to_length, contig_to_kegg_counter, contigs_fasta, diamond_result_file)
    
    ### Use contig index instead of contig name to save memory
    contig_to_index, index_to_contig = contig_manager.make_contig_index(contigs_in_bins)

    contig_to_kegg_counter = contig_manager.apply_contig_index(contig_to_index, contig_to_kegg_counter)
    contig_to_genes = contig_manager.apply_contig_index(contig_to_index, contig_to_genes)
    contig_to_length = contig_manager.apply_contig_index(contig_to_index, contig_to_length) 
    
    bin_manager.rename_bin_contigs(original_bins, contig_to_index)

    ## Extract cds metadata ##

    logging.info('Compute cds metadata.')
    contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length = cds.get_contig_cds_metadata(contig_to_genes, threads) 
    
    contig_info = {'contig_to_cds_count':contig_to_cds_count,
                    'contig_to_aa_counter':contig_to_aa_counter,
                    'contig_to_aa_length':contig_to_aa_length,
                    'contig_to_kegg_counter':contig_to_kegg_counter,
                    'contig_to_length':contig_to_length}

    logging.info('Add size and assess quality of input bins')

    # TODO paralellize
    # original_bins = bin_quality.add_bin_metrics_in_parallel(original_bins, contig_info, threads)

    bin_quality.add_bin_metrics(original_bins, contig_info, contamination_weigth)
        
    logging.info('Create intermediate bins:')
    new_bins = bin_manager.create_intermediate_bins(bin_set_name_to_bins)
 
    logging.info(f'Assess quality for supplementary intermediate bins.')
    new_bins = bin_quality.add_bin_metrics(new_bins, contig_info, threads, contamination_weigth)

    
    # bin_quality.add_bin_size_and_N50(new_bins, contig_to_length)

    # postProcessor = modelPostprocessing.modelProcessor(1)
    # logging.info(f'Assess bin quality of {len(new_bins)} new bins created from intersection, difference or unions in bin graph.')
    # bin_quality.assess_bins_quality_by_chunk(new_bins, contig_to_kegg_counter, contig_to_cds_count,
    #                                 contig_to_aa_counter, contig_to_aa_length,
    #                                 postProcessor=postProcessor,  threads=threads)
        
   
    logging.info('Dereplicating input bins and new bins')
    all_bins = original_bins | new_bins

    # if debug:
    #     import pickle
    #     # import networkx as nx
    #     with open(os.path.join(out_tmp_dir, 'all_bin.p' ), "bw") as fl:
    #         pickle.dump( all_bins, fl)
    #     logging.debug('Writting all bins info ')
    #     # write_bin_info_debug(all_bins, os.path.join(out_tmp_dir, 'all_bins.tsv'))

    #     # G = bin_manager.get_bin_graph_with_attributes(all_bins, contig_to_length)

    #     # nx.write_edgelist(G, os.path.join(out_tmp_dir, "bin_graph_edglist"))
    #     # # df = nx.to_pandas_adjacency(G)
    #     # # df.to_csv(os.path.join(out_tmp_dir, "bin_graph_edglist.tsv"), sep='\t')
    #     # nx.write_weighted_edgelist(G, os.path.join(out_tmp_dir, "bin_graph_w_edglist.tsv"), delimiter='\t')

    #     # with open(os.path.join(out_tmp_dir, "bin_graph_edglist"), "w") as fl:
    #     #     for e in nx.edges(G):

    #     #         print(e)
    #     #         print(dir(e))
    #     #         print(e)


    logging.info('Select best bins')
    selected_bins = bin_manager.select_best_bins(all_bins)

    logging.info(f'Filtering bins: only bins with completeness >= {min_completeness} are kept')
    selected_bins = [b for b in selected_bins if b.completeness >= min_completeness]

    logging.info(f'Writing selected bins in {final_bin_report}')

    for b in selected_bins:
        b.contigs = {index_to_contig[c_index] for c_index in b.contigs }
    write_bin_info(selected_bins, final_bin_report)

    write_bins_fasta(selected_bins, contigs_fasta, outdir_final_bin_set)

    if debug:
        for sb in selected_bins:
            if sb.completeness >= hq_min_completeness and sb.contamination <= hq_max_conta:
                    logging.debug(f"{sb}, {sb.completeness}, {sb.contamination}")


    hq_bins = len([sb for sb in selected_bins if sb.completeness >= hq_min_completeness and sb.contamination <= hq_max_conta ]) 
    hq_bins_single = len([sb for sb in selected_bins if sb.completeness >= hq_min_completeness and sb.contamination <=  hq_max_conta  and len(sb.contigs) ==1 ]) 
    logging.info(f'{hq_bins}/{len(selected_bins)} selected bins have a high quality (completeness >= {hq_min_completeness} and contamination <= {hq_max_conta}).')
    logging.info(f'{hq_bins_single}/{len(selected_bins)} selected bins have a high quality and are made of only one contig.')

# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
