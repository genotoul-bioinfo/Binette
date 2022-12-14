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
import bin_manager
import file_manager
import os
import cds
import diamond
import bin_quality
from checkm2 import modelPostprocessing

# import pkg_resources


# try:
#     PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
# except pkg_resources.DistributionNotFound:
#     PROGRAM_VERSION = "undefined_version"


def init_logging(verbose):
    '''Initialise logging.'''

    if verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt="[%Y-%m-%d %H:%M:%S]")

    else:
        logging.basicConfig(format='%(asctime)s %(levelname)s - %(message)s',
                    datefmt="[%Y-%m-%d %H:%M:%S]")

    logging.info('Program started')
    logging.info(f'command line: {" ".join(sys.argv)}', )



def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(description="...",
                            formatter_class=ArgumentDefaultsHelpFormatter)
                        
    parser.add_argument("-i", "--bins", nargs='+', required=True, 
                        help="Bin folders containing each bin in a fasta file.")

    parser.add_argument("-c", "--contigs", required=True, help="Contigs in fasta format.")
    parser.add_argument("-t", "--threads", default=1, type=int, help="Number of threads.")

    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")

    # parser.add_argument('--version',
    #                     action='version',
    #                     version='%(prog)s ' + PROGRAM_VERSION)

    args = parser.parse_args()
    return args

def infer_bin_name_from_bin_dir(bin_dirs):
    # remove common prefix of bin dir to get nicer label
    commonprefix_len = len(os.path.commonprefix(bin_dirs))
    bin_name_to_bin_dir = {d[commonprefix_len:]:d for d in bin_dirs}
    return bin_name_to_bin_dir

def write_bin_info(bins, output):

    header = ['origin', "name", 'completeness', 'contamination', "score", 'size', 'contig_count',  'contigs']
    with open(output, "w") as fl:
        fl.write('\t'.join(header)+"\n")
        for bin_obj in bins:

            line = [bin_obj.origin, bin_obj.name,
                    bin_obj.completeness, bin_obj.contamination, bin_obj.score, 
                     len(bin_obj.contigs), bin_obj.length,
                     ";".join(bin_obj.contigs)]

            fl.write("\t".join((str(e) for e in line)) + '\n')

def check_contig_consistency(contigs_from_assembly, contigs_from_elsewhere, assembly_file, elsewhere_file):
    are_contigs_consistent = len(set(contigs_from_elsewhere) | set(contigs_from_assembly)) <= len(set(contigs_from_assembly))
    message = f"{len(set(contigs_from_elsewhere) - set(contigs_from_assembly))} contigs found in file {elsewhere_file} were not found in assembly_file ({assembly_file})."
    assert are_contigs_consistent, message

def main():
    "Orchestrate the execution of the program"

    args = parse_arguments()

    init_logging(args.verbose)

    bin_dirs = args.bins
    contigs_fasta = args.contigs
    threads = args.threads
    faa_file = 'tmp_head.faa'
    diamond_result_file = 'diamond_result.tsv' 
    run_tool = False


    logging.info('Parse bin directories.')
    bin_name_to_bin_dir = infer_bin_name_from_bin_dir(bin_dirs)
    bin_name_to_bins = bin_manager.parse_bin_directories(bin_name_to_bin_dir)
    original_bins = bin_manager.dereplicate_bin_sets(bin_name_to_bins.values()) 
    contigs_in_bins = bin_manager.get_contigs_in_bins(original_bins) 

    logging.info('Parse fasta file of assembly.')
    contigs_object = file_manager.parse_fasta_file(contigs_fasta)
    contig_to_length = {seq.name:len(seq) for seq in contigs_object if contigs_in_bins}


    if run_tool:
        contig_to_genes = cds.predict(contigs_fasta, faa_file, threads)
    else:
        logging.info('Parse faa file.')
        contig_to_genes = cds.parse_faa_file(faa_file)
        check_contig_consistency(contig_to_length, contig_to_genes, contigs_fasta, faa_file)

    # from genes extract contig cds metadata
    logging.info('Compute cds metadata.')
    contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length = cds.get_contig_cds_metadata(contig_to_genes) 

    if run_tool:
        diamond_db_path = diamond.get_checkm2_db()
        diamond.run(faa_file, diamond_result_file, diamond_db_path, threads)

    contig_to_kegg_counter = bin_quality.get_contig_to_kegg_id(diamond_result_file)

    ## Check contigs from diamond vs input assembly consistency
    check_contig_consistency(contig_to_length, contig_to_kegg_counter, contigs_fasta, diamond_result_file)
    

    logging.info('Add size and assess quality of input bins')
    postProcessor = modelPostprocessing.modelProcessor(threads)

    for bin_set_id, bins in bin_name_to_bins.items():
        logging.info(f'{bin_set_id} - {len(bins)} ')
        logging.info('Add size')
        bin_manager.add_bin_size(bins, contig_to_length)
        logging.info('Asses quality')
        bin_quality.assess_bins_quality(bins, contig_to_kegg_counter, contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length, postProcessor=postProcessor,  threads=threads)
        

    logging.info('Create intermediate bins:')

    logging.info('Making bin graph...')
    connected_bins_graph = bin_manager.from_bin_sets_to_bin_graph(bin_name_to_bins)

    logging.info('Create intersection bins...')
    intersection_bins = bin_manager.get_intersection_bins(connected_bins_graph)
    logging.info(f'{len(intersection_bins)} bins created on intersections.')

    logging.info('Create difference bin...')
    difference_bins =  bin_manager.get_difference_bins(connected_bins_graph)
    logging.info(f'{len(difference_bins)} bins created based on symetric difference.')


    logging.info('Create get_union_bins bin...')
    union_bins =  bin_manager.get_union_bins(connected_bins_graph)
    logging.info(f'{len(union_bins)} bins created on unions.')

    new_bins = difference_bins | intersection_bins | union_bins 

    logging.info('Assess bin quality of new bins created from intersection, difference or unions in bin graph.')
    bin_quality.assess_bins_quality(new_bins, contig_to_kegg_counter, contig_to_cds_count,
                                    contig_to_aa_counter, contig_to_aa_length,
                                    postProcessor=postProcessor,  threads=threads)
        

    logging.info('Dereplicating input bins and new bins')
    bin_sets = bin_name_to_bins.values()
    original_bins = bin_manager.dereplicate_bin_sets(bin_sets)
    all_bins = original_bins | new_bins


    logging.info('Select best bin')
    selected_bins = bin_manager.select_best_bins(all_bins)

    logging.info('Writing selected bins')
    write_bin_info(selected_bins, 'selected_bins.tsv')

    # import pickle
    # with open("all_bin.p", "bw") as fl:
    #     pickle.dump( all_bins, fl)





# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
