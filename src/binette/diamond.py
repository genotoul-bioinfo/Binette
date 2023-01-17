import subprocess
import logging
import sys
import shutil
import re
import sys
import pandas as pd
from collections import Counter

from checkm2 import keggData


def get_checkm2_db():

    if shutil.which('checkm2') is None:
        logging.error("Make sure checkm2 is on your system path.")
        sys.exit(1)

    checkm2_database_raw = subprocess.run(['checkm2', 'database', '--current'], text=True,  stderr=subprocess.PIPE)

    if checkm2_database_raw.returncode != 0:
        logging.error(f'Something went wrong with checkm2:\n=======\n{checkm2_database_raw.stderr}========')
        sys.exit(1)


    reg_result = re.search("INFO: (/.*.dmnd)", checkm2_database_raw.stderr)
    
    try:
        db_path = reg_result.group(1)

    except AttributeError:
        logging.error(f'Something went wrong when retrieving checkm2 db path:\n{checkm2_database_raw.stderr}')
        sys.exit(1)

    return db_path


def check_diamond_exists():
    """Check to see if Diamond is on the system before we try to run it."""

    # Assume that a successful diamond help returns 0 and anything
    # else returns something non-zero

    if shutil.which('diamond') is None:
        logging.error("Make sure diamond is on your system path.")
        sys.exit(1)
    

def run(faa_file, output, db, log, threads=1, query_cover=80, subject_cover=80, percent_id=30, evalue=1e-05, low_mem=False):
    
    check_diamond_exists()


    blocksize = 0.5 if low_mem else 2

    cmd = "diamond blastp --outfmt 6 --max-target-seqs 1 " \
        f"--query {faa_file} " \
        f"-o {output} " \
        f"--threads {threads} " \
        f"--db {db} " \
        f"--query-cover {query_cover} " \
        f"--subject-cover {subject_cover} " \
        f"--id {percent_id} " \
        f"--evalue {evalue} --block-size {blocksize} 2> {log}"


    logging.info('Running diamond')
    logging.info(cmd)

    
    run = subprocess.run(cmd, shell=True)

    if run.returncode != 0:
        logging.error(f'An error occured while running DIAMOND. check log file: {log}')
        sys.exit(1)

    logging.info('Finished Running DIAMOND')


def get_contig_to_kegg_id(diamond_result_file):

    diamon_results_df = pd.read_csv(diamond_result_file, sep='\t', usecols=[0, 1], names=['ProteinID', 'annotation'])
    diamon_results_df[['Ref100_hit', 'Kegg_annotation']] = diamon_results_df['annotation'].str.split('~', n=1, expand=True)
    diamon_results_df

    ''' Get a list of default KO id's from data
        Available categories are the keys in DefaultValues.feature_ordering
        Here, returns an ordered set of KEGG ID's and sets to 0 
    '''
    KeggCalc = keggData.KeggCalculator()
    defaultKOs = KeggCalc.return_default_values_from_category('KO_Genes')

    #Remove from diamon_results_df any KOs not currently used by checkm2
    diamon_results_df = diamon_results_df.loc[diamon_results_df['Kegg_annotation'].isin(defaultKOs.keys())]
    diamon_results_df['contig'] = diamon_results_df['ProteinID'].str.split('_', n=-1).str[:-1].str.join('_')
    #diamon_results_df[diamon_results_df['Kegg_annotation']]
    # group by contig and create a counter with kegg_annotation
    contig_to_kegg_counter = diamon_results_df.groupby("contig").agg({'Kegg_annotation':Counter}).reset_index() # ['Kegg_annotation'].apply(Counter)

    # create a simple dict with contig --> kegg_counter
    contig_to_kegg_counter = dict(zip(contig_to_kegg_counter['contig'], contig_to_kegg_counter['Kegg_annotation']))

    return contig_to_kegg_counter