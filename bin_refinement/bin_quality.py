#!/usr/bin/env python3

"""

"""


import os 
from collections import Counter
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

@measure_memory
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

@measure_memory
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

@measure_memory
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
    


@measure_memory
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

