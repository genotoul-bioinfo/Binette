#!/usr/bin/env python3
import logging
import os
from collections import Counter
from itertools import islice
from typing import Dict, Iterable, Optional, Tuple, Iterator, Set

import numpy as np
import pandas as pd
from binette.bin_manager import Bin


# Suppress unnecessary TensorFlow warnings
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
logging.getLogger("tensorflow").setLevel(logging.FATAL)


from checkm2 import keggData, modelPostprocessing, modelProcessing  # noqa: E402


def get_bins_metadata_df(bins: Iterable[Bin], contig_to_cds_count: Dict[str, int], contig_to_aa_counter: Dict[str, Counter], contig_to_aa_length: Dict[str, int]) -> pd.DataFrame:
    """
    Generate a DataFrame containing metadata for a list of bins.

    :param bins: A list of bin objects.
    :param contig_to_cds_count: A dictionary mapping contig names to CDS counts.
    :param contig_to_aa_counter: A dictionary mapping contig names to amino acid composition (Counter object).
    :param contig_to_aa_length: A dictionary mapping contig names to total amino acid length.
    :return: A DataFrame containing bin metadata.
    """

    metadata_order = keggData.KeggCalculator().return_proper_order("Metadata")

    # Get bin metadata
    bin_metadata_list = []
    for bin_obj in bins:
        bin_metadata = {
            "Name": bin_obj.id,
            "CDS": sum((contig_to_cds_count[c] for c in bin_obj.contigs if c in contig_to_cds_count)),
            "AALength": sum((contig_to_aa_length[c] for c in bin_obj.contigs if c in contig_to_aa_length)),
        }

        bin_aa_counter = Counter()
        for contig in bin_obj.contigs:
            if contig in contig_to_aa_counter:
                bin_aa_counter += contig_to_aa_counter[contig]

        bin_metadata.update(dict(bin_aa_counter))
        bin_metadata_list.append(bin_metadata)

    metadata_df = pd.DataFrame(bin_metadata_list, columns=["Name"] + metadata_order)

    metadata_df = metadata_df.fillna(0)

    metadata_df = metadata_df.astype({col: int for col in metadata_order})

    metadata_df = metadata_df.set_index("Name", drop=False)
    return metadata_df

def get_diamond_feature_per_bin_df(bins: Iterable[Bin], contig_to_kegg_counter: Dict[str, Counter]) -> Tuple[pd.DataFrame, int]:
    """
    Generate a DataFrame containing Diamond feature counts per bin and completeness information for pathways, categories, and modules.

    :param bins: A list of bin objects.
    :param contig_to_kegg_counter: A dictionary mapping contig names to KEGG annotation counters.
    :type bins: List
    :type contig_to_kegg_counter: Dict[str, Counter]
    :return: A tuple containing the DataFrame and the number of default KEGG orthologs (KOs).
    :rtype: Tuple[pd.DataFrame, int]
    """
    KeggCalc = keggData.KeggCalculator()
    defaultKOs = KeggCalc.return_default_values_from_category("KO_Genes")

    bin_to_ko_counter = {}
    for bin_obj in bins:
        bin_ko_counter = Counter()
        for contig in bin_obj.contigs:
            try:
                bin_ko_counter += contig_to_kegg_counter[contig]
            except KeyError:
                # No KO annotation found in this contig
                continue

        bin_to_ko_counter[bin_obj.id] = bin_ko_counter

    ko_count_per_bin_df = pd.DataFrame(bin_to_ko_counter, index=defaultKOs).transpose().fillna(0)
    ko_count_per_bin_df = ko_count_per_bin_df.astype(int)
    ko_count_per_bin_df["Name"] = ko_count_per_bin_df.index

    logging.info("Calculating completeness of pathways and modules.")
    logging.debug("Calculating pathway completeness information")
    KO_pathways = KeggCalc.calculate_KO_group("KO_Pathways", ko_count_per_bin_df.copy())

    logging.debug("Calculating category completeness information")
    KO_categories = KeggCalc.calculate_KO_group("KO_Categories", ko_count_per_bin_df.copy())

    logging.debug("Calculating module completeness information")
    KO_modules = KeggCalc.calculate_module_completeness(ko_count_per_bin_df.copy())

    diamond_complete_results = pd.concat([ko_count_per_bin_df, KO_pathways, KO_modules, KO_categories], axis=1)

    return diamond_complete_results, len(defaultKOs)


def compute_N50(list_of_lengths) -> int:
    """
    Calculate N50 for a sequence of numbers.

    :param list_of_lengths: List of numbers.
    :param list_of_lengths: list
    :return: N50 value.
    """
    list_of_lengths = sorted(list_of_lengths)
    sum_len = sum(list_of_lengths)

    cum_length = 0
    length = 0
    for length in list_of_lengths:
        if cum_length + length >= sum_len / 2:
            return length
        cum_length += length
    return length

def add_bin_size_and_N50(bins: Iterable[Bin], contig_to_size: Dict[str,int]):
    """
    Add bin size and N50 to a list of bin objects.

    :param bins: List of bin objects.
    :param contig_to_size: Dictionary mapping contig names to their sizes.
    """
    for bin_obj in bins:
        lengths = [contig_to_size[c] for c in bin_obj.contigs]
        n50 = compute_N50(lengths)

        bin_obj.add_length(sum(lengths))
        bin_obj.add_N50(n50)


def add_bin_metrics(bins: Set[Bin], contig_info: Dict, contamination_weight: float, threads: int = 1):
    """
    Add metrics to a Set of bins.

    :param bins: Set of bin objects.
    :param contig_info: Dictionary containing contig information.
    :param contamination_weight: Weight for contamination assessment.
    :param threads: Number of threads for parallel processing (default is 1).

    :return: List of processed bin objects.
    """
    postProcessor = modelPostprocessing.modelProcessor(threads)

    contig_to_kegg_counter = contig_info["contig_to_kegg_counter"]
    contig_to_cds_count = contig_info["contig_to_cds_count"]
    contig_to_aa_counter = contig_info["contig_to_aa_counter"]
    contig_to_aa_length = contig_info["contig_to_aa_length"]
    contig_to_length = contig_info["contig_to_length"]

    logging.info("Assess bin length and N50")

    add_bin_size_and_N50(bins, contig_to_length)

    logging.info("Assess bin quality")
    assess_bins_quality_by_chunk(
        bins,
        contig_to_kegg_counter,
        contig_to_cds_count,
        contig_to_aa_counter,
        contig_to_aa_length,
        contamination_weight,
        postProcessor,
    )
    return bins


def chunks(iterable: Iterable, size: int) -> Iterator[Tuple]:
    """
    Generate adjacent chunks of data from an iterable.

    :param iterable: The iterable to be divided into chunks.
    :param size: The size of each chunk.
    :return: An iterator that produces tuples of elements in chunks.
    """
    it = iter(iterable)
    return iter(lambda: tuple(islice(it, size)), ())


def assess_bins_quality_by_chunk(bins: Iterable[Bin],
    contig_to_kegg_counter: Dict,
    contig_to_cds_count: Dict,
    contig_to_aa_counter: Dict,
    contig_to_aa_length: Dict,
    contamination_weight: float,
    postProcessor:Optional[modelPostprocessing.modelProcessor] = None,
    threads: int = 1,
    chunk_size: int = 2500):
    """
    Assess the quality of bins in chunks.

    This function assesses the quality of bins in chunks to improve processing efficiency.

    :param bins: List of bin objects.
    :param contig_to_kegg_counter: Dictionary mapping contig names to KEGG counters.
    :param contig_to_cds_count: Dictionary mapping contig names to CDS counts.
    :param contig_to_aa_counter: Dictionary mapping contig names to amino acid counters.
    :param contig_to_aa_length: Dictionary mapping contig names to amino acid lengths.
    :param contamination_weight: Weight for contamination assessment.
    :param postProcessor: post-processor from checkm2
    :param threads: Number of threads for parallel processing (default is 1).
    :param chunk_size: The size of each chunk. 
    """

    for i, chunk_bins_iter in enumerate(chunks(bins, chunk_size)):
        chunk_bins = set(chunk_bins_iter)
        logging.debug(f"chunk {i}: assessing quality of {len(chunk_bins)}")
        assess_bins_quality(
            bins=chunk_bins,
            contig_to_kegg_counter= contig_to_kegg_counter,
            contig_to_cds_count=contig_to_cds_count,
            contig_to_aa_counter=contig_to_aa_counter,
            contig_to_aa_length=contig_to_aa_length,
            contamination_weight=contamination_weight,
            postProcessor=postProcessor,
            threads=threads
        )

def assess_bins_quality(
    bins: Iterable[Bin],
    contig_to_kegg_counter: Dict,
    contig_to_cds_count: Dict,
    contig_to_aa_counter: Dict,
    contig_to_aa_length: Dict,
    contamination_weight: float,
    postProcessor: Optional[modelPostprocessing.modelProcessor] = None,
    threads: int = 1,):
    """
    Assess the quality of bins.

    This function assesses the quality of bins based on various criteria and assigns completeness and contamination scores. 
    This code is taken from checkm2 and adjusted

    :param bins: List of bin objects.
    :param contig_to_kegg_counter: Dictionary mapping contig names to KEGG counters.
    :param contig_to_cds_count: Dictionary mapping contig names to CDS counts.
    :param contig_to_aa_counter: Dictionary mapping contig names to amino acid counters.
    :param contig_to_aa_length: Dictionary mapping contig names to amino acid lengths.
    :param contamination_weight: Weight for contamination assessment.
    :param postProcessor: A post-processor from checkm2
    :param threads: Number of threads for parallel processing (default is 1).
    """
    if postProcessor is None:
        postProcessor = modelPostprocessing.modelProcessor(threads)

    metadata_df = get_bins_metadata_df(bins, contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length)

    diamond_complete_results, ko_list_length = get_diamond_feature_per_bin_df(bins, contig_to_kegg_counter)
    diamond_complete_results = diamond_complete_results.drop(columns=["Name"])

    feature_vectors = pd.concat([metadata_df, diamond_complete_results], axis=1)
    feature_vectors = feature_vectors.sort_values(by="Name")

    # 4: Call general model & specific models and derive predictions"""
    modelProc = modelProcessing.modelProcessor(threads)

    vector_array = feature_vectors.iloc[:, 1:].values.astype(np.float)

    logging.info("Predicting completeness and contamination using the general model.")
    general_results_comp, general_results_cont = modelProc.run_prediction_general(vector_array)

    logging.info("Predicting completeness using the specific model.")
    specific_model_vector_len = (
        ko_list_length + len(metadata_df.columns)
    ) - 1

    # also retrieve scaled data for CSM calculations
    specific_results_comp, scaled_features = modelProc.run_prediction_specific(vector_array, specific_model_vector_len)

    logging.info("Using cosine similarity to reference data to select an appropriate predictor model.")

    final_comp, final_cont, models_chosen, csm_array = postProcessor.calculate_general_specific_ratio(
        vector_array[:, 20], scaled_features, general_results_comp, general_results_cont, specific_results_comp
    )

    final_results = feature_vectors[["Name"]].copy()
    final_results["Completeness"] = np.round(final_comp, 2)
    final_results["Contamination"] = np.round(final_cont, 2)

    for bin_obj in bins:
        completeness = final_results.at[bin_obj.id, "Completeness"]
        contamination = final_results.at[bin_obj.id, "Contamination"]

        bin_obj.add_quality(completeness, contamination, contamination_weight)
