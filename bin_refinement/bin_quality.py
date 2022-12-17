#!/usr/bin/env python3

"""

"""


import os 
from collections import Counter
from itertools import islice
import pandas as pd
import numpy as np
import logging


# For unnessesary tensorflow warnings:
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
logging.getLogger('tensorflow').setLevel(logging.FATAL)

from checkm2 import keggData
from checkm2 import modelProcessing
from checkm2 import modelPostprocessing


from memory_control import measure_memory


def get_bins_metadata_df(bins, contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length):

    metadata_order = keggData.KeggCalculator().return_proper_order('Metadata')
    ## Get bin metadata
    bin_metadata_list = []
    for bin_obj in bins:
        bin_metadata = {"Name":bin_obj.id,
                        'CDS': sum((contig_to_cds_count[c] for c in bin_obj.contigs if c in contig_to_cds_count)), 
                        "AALength":sum((contig_to_aa_length[c] for c in bin_obj.contigs if c in contig_to_aa_length)),
                        }

        bin_aa_counter  = Counter()
        for contig in bin_obj.contigs:
            if contig in contig_to_aa_counter:
                bin_aa_counter += contig_to_aa_counter[contig]

        bin_metadata.update(dict(bin_aa_counter))

        bin_metadata_list.append(bin_metadata)

    metadata_order = keggData.KeggCalculator().return_proper_order('Metadata')

    metadata_df = pd.DataFrame(bin_metadata_list, columns = ["Name"] + metadata_order)

    metadata_df = metadata_df.fillna(0)

    metadata_df = metadata_df.astype({col: int for col in metadata_order})
    
    metadata_df = metadata_df.set_index('Name', drop=False)
    return metadata_df


def get_diamond_feature_per_bin_df(bins, contig_to_kegg_counter):

    KeggCalc = keggData.KeggCalculator()
    defaultKOs = KeggCalc.return_default_values_from_category('KO_Genes')

    bin_to_ko_counter = {}
    for bin_obj in bins:
        bin_ko_counter = Counter()
        for contig in bin_obj.contigs:
            try:
                bin_ko_counter += contig_to_kegg_counter[contig]
            except KeyError:
                # no ko annotation found in this contig
                continue

        bin_to_ko_counter[bin_obj.id] = bin_ko_counter

    ko_count_per_bin_df = pd.DataFrame(bin_to_ko_counter, index=defaultKOs).transpose().fillna(0)
    ko_count_per_bin_df = ko_count_per_bin_df.astype(int)
    ko_count_per_bin_df['Name'] = ko_count_per_bin_df.index


    logging.info('Calculating completeness of pathways and modules.')
    logging.debug('Calculating pathway completeness information')
    KO_pathways = KeggCalc.calculate_KO_group('KO_Pathways', ko_count_per_bin_df.copy())


    logging.debug('Calculating category completeness information')
    KO_categories = KeggCalc.calculate_KO_group('KO_Categories', ko_count_per_bin_df.copy())

    logging.debug('Calculating module completeness information')
    KO_modules = KeggCalc.calculate_module_completeness(ko_count_per_bin_df.copy())

    diamond_complete_results = pd.concat([ko_count_per_bin_df, KO_pathways, KO_modules, KO_categories], axis=1)

    return diamond_complete_results, len(defaultKOs)
    

def chunks(iterable, size):
    """Generate adjacent chunks of data"""
    it = iter(iterable)
    return iter(lambda: tuple(islice(it, size)), ())

def assess_bins_quality_by_chunk(bins, contig_to_kegg_counter, contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length, postProcessor=None, threads=1):
    n = 2500
    
    for i, chunk_bins_iter in enumerate(chunks(bins, n)):
        chunk_bins = set(chunk_bins_iter)
        logging.debug(f'chunk {i}: assessing quality of {len(chunk_bins)}')
        assess_bins_quality(chunk_bins, contig_to_kegg_counter, contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length, postProcessor)


def assess_bins_quality(bins, contig_to_kegg_counter, contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length, postProcessor=None, threads=1):

    if postProcessor is None:
        postProcessor = modelPostprocessing.modelProcessor(threads) 


    metadata_df = get_bins_metadata_df(bins, contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length)
    
    diamond_complete_results, ko_list_length = get_diamond_feature_per_bin_df(bins, contig_to_kegg_counter)
    diamond_complete_results = diamond_complete_results.drop(columns=['Name'])

    feature_vectors = pd.concat([metadata_df, diamond_complete_results], axis=1)
    feature_vectors = feature_vectors.sort_values(by='Name')
    

    ''' 4: Call general model & specific models and derive predictions'''
    modelProc = modelProcessing.modelProcessor(threads)

    vector_array = feature_vectors.iloc[:, 1:].values.astype(np.float)

    logging.info('Predicting completeness and contamination using general model.')
    general_results_comp, general_results_cont = modelProc.run_prediction_general(vector_array)


    logging.info('Predicting completeness using specific model.')
    specific_model_vector_len = (ko_list_length + len(metadata_df.columns)) - 1  # -1 = without name TODO a bit ugly - maybe just calculate length on setup somewhere

    # also retrieve scaled data for CSM calculations
    specific_results_comp, scaled_features = modelProc.run_prediction_specific(vector_array, specific_model_vector_len)

    logging.info('Using cosine simlarity to reference data to select appropriate predictor model.')

    # postProcessor = modelPostprocessing.modelProcessor(threads)
    final_comp, final_cont, models_chosen, csm_array = postProcessor.calculate_general_specific_ratio(
        vector_array[:, 20], 
        scaled_features,
        general_results_comp,
        general_results_cont,
        specific_results_comp)


    final_results = feature_vectors[['Name']].copy()
    final_results['Completeness'] = np.round(final_comp, 2)
    final_results['Contamination'] = np.round(final_cont, 2)

    for bin_obj in bins:
        completeness = final_results.loc[bin_obj.id, 'Completeness']
        contamination = final_results.loc[bin_obj.id, 'Contamination']

        bin_obj.add_quality(completeness, contamination)

