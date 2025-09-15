#!/usr/bin/env python3
import logging
import os
import sys
from collections import Counter, defaultdict
from itertools import islice
from typing import Dict, Iterable, Tuple, Iterator, List

import lazy_loader as lazy
import numpy as np
import pandas as pd
from tqdm import tqdm

from binette.bin_manager import Bin
from checkm2 import keggData

# Suppress unnecessary TensorFlow warnings
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
logging.getLogger("tensorflow").setLevel(logging.FATAL)

# Lazy loaders for checkm2 components that import keras
# These will only be imported when explicitly called
_modelPostprocessing = None
_modelProcessing = None


def get_modelPostprocessing():
    """Lazy load modelPostprocessing module only when needed"""
    global _modelPostprocessing
    if _modelPostprocessing is None:
        # Only import keras when absolutely needed
        sys.modules["keras"] = lazy.attach("keras")
        sys.modules["keras.models"] = lazy.attach("keras.models")
        from checkm2 import modelPostprocessing

        _modelPostprocessing = modelPostprocessing
    return _modelPostprocessing


def get_modelProcessing():
    """Lazy load modelProcessing module only when needed"""
    global _modelProcessing
    if _modelProcessing is None:
        # Only import keras when absolutely needed
        sys.modules["keras"] = lazy.attach("keras")
        sys.modules["keras.models"] = lazy.attach("keras.models")
        from checkm2 import modelProcessing

        _modelProcessing = modelProcessing
    return _modelProcessing


def get_bins_metadata_df(
    bins: List[Bin],
    contig_to_cds_count: Dict[str, int],
    contig_to_aa_counter: Dict[str, Counter],
    contig_to_aa_length: Dict[str, int],
) -> pd.DataFrame:
    """
    Optimized: Generate a DataFrame containing metadata for a list of bins.
    Handles contigs that appear in multiple bins.
    """
    metadata_order = keggData.KeggCalculator().return_proper_order("Metadata")
    bin_keys = [b.contigs_key for b in bins]

    # --- Pre-aggregate CDS and AA length ---
    cds_per_bin = defaultdict(int)
    aa_len_per_bin = defaultdict(int)
    aa_counter_per_bin = defaultdict(Counter)

    # map contigs → all bins they belong to
    contig_to_bins = defaultdict(list)
    for b in bins:
        for c in b.contigs:
            contig_to_bins[c].append(b.contigs_key)

    # distribute CDS counts
    for contig, cds in contig_to_cds_count.items():
        for bin_key in contig_to_bins.get(contig, []):
            cds_per_bin[bin_key] += cds

    # distribute AA lengths
    for contig, length in contig_to_aa_length.items():
        for bin_key in contig_to_bins.get(contig, []):
            aa_len_per_bin[bin_key] += length

    # distribute AA counters
    for contig, counter in contig_to_aa_counter.items():
        for bin_key in contig_to_bins.get(contig, []):
            aa_counter_per_bin[bin_key].update(counter)

    # --- Build rows ---
    rows = []
    for key in bin_keys:
        row = {
            "Name": key,
            "CDS": cds_per_bin.get(key, 0),
            "AALength": aa_len_per_bin.get(key, 0),
        }
        row.update(aa_counter_per_bin.get(key, {}))
        rows.append(row)

    # --- Construct DataFrame directly ---
    metadata_df = pd.DataFrame(rows).fillna(0)

    # Ensure column order
    all_cols = ["Name"] + metadata_order
    for col in metadata_order:
        if col not in metadata_df.columns:
            metadata_df[col] = 0

    metadata_df = metadata_df[all_cols].astype({col: int for col in metadata_order})
    metadata_df = metadata_df.set_index("Name", drop=False)

    return metadata_df


def get_diamond_feature_per_bin_df(
    bins: List[Bin], contig_to_kegg_counter: Dict[str, Counter]
) -> Tuple[pd.DataFrame, int]:
    """
    Optimized: Generate a DataFrame containing Diamond feature counts per bin,
    including KEGG KO counts and completeness information for pathways, categories, and modules.
    Handles contigs that may belong to multiple bins.
    """
    KeggCalc = keggData.KeggCalculator()
    defaultKOs = KeggCalc.return_default_values_from_category("KO_Genes")
    bin_keys = [b.contigs_key for b in bins]

    # --- Build contig → bins mapping ---
    contig_to_bins = defaultdict(list)
    for b in bins:
        for c in b.contigs:
            contig_to_bins[c].append(b.contigs_key)

    # --- Aggregate KO counters per bin ---
    bin_to_ko_counter = {}
    for bin_obj in bins:
        bin_ko_counter = Counter()
        for contig in bin_obj.contigs:
            ko_counter = contig_to_kegg_counter.get(contig)
            if ko_counter:
                bin_ko_counter.update(ko_counter)
        bin_to_ko_counter[bin_obj.contigs_key] = bin_ko_counter

    # --- Build KO count DataFrame directly ---
    ko_count_per_bin_df = (
        pd.DataFrame.from_dict(bin_to_ko_counter, orient="index")
        .reindex(bin_keys)  # keep bin order
        .fillna(0)
        .astype(int)
    )

    # Ensure all defaultKOs exist
    ko_count_per_bin_df = ko_count_per_bin_df.reindex(
        columns=list(defaultKOs), fill_value=0
    )

    # ko_count_per_bin_df.index.name = "Name"
    ko_count_per_bin_df["Name"] = ko_count_per_bin_df.index

    # --- Calculate higher-level completeness ---
    logging.debug("Calculating completeness of pathways, categories, and modules.")
    KO_pathways = calculate_KO_group(KeggCalc, "KO_Pathways", ko_count_per_bin_df)
    KO_categories = calculate_KO_group(KeggCalc, "KO_Categories", ko_count_per_bin_df)

    KO_modules = calculate_module_completeness(KeggCalc, ko_count_per_bin_df)

    # --- Concatenate results ---
    diamond_complete_results = pd.concat(
        [ko_count_per_bin_df, KO_pathways, KO_modules, KO_categories], axis=1
    )

    return diamond_complete_results, len(defaultKOs)


def calculate_KO_group(
    KeggCalc: keggData.KeggCalculator, group: str, KO_gene_data: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate the completeness of KEGG feature groups per bin.

    :param KeggCalc: An instance of KeggCalculator containing KEGG mappings.
    :param group: Feature group name (e.g., "KO_Pathways", "KO_Categories").
    :param KO_gene_data: DataFrame containing KO counts per bin with last column "Name".

    :return: DataFrame with completeness values for each feature vector in the group.
    """

    # last column is 'Name'
    data = KO_gene_data.drop(columns=["Name"]).values
    n_bins = data.shape[0]

    # Build output DataFrame
    ordered_entries = KeggCalc.return_default_values_from_category(group)
    feature_vectors = list(ordered_entries.keys())
    n_features = len(feature_vectors)

    # Create empty numpy array for results
    result = np.zeros((n_bins, n_features), dtype=float)

    # Map Kegg_IDs to column indices in KO_gene_data
    col_map = {ko: idx for idx, ko in enumerate(KO_gene_data.columns[:-1])}

    for f_idx, vector in enumerate(feature_vectors):
        # KOs belonging to this feature vector
        kegg_ids = KeggCalc.path_category_mapping.loc[
            KeggCalc.path_category_mapping[group] == vector, "Kegg_ID"
        ].values

        # Only keep KOs present in DataFrame columns
        present_cols = [col_map[ko] for ko in kegg_ids if ko in col_map]
        if not present_cols:
            continue

        # Presence/absence: values >1 -> 1
        vals = data[:, present_cols]
        vals[vals > 1] = 1
        result[:, f_idx] = vals.sum(axis=1) / len(kegg_ids)

    return pd.DataFrame(result, columns=feature_vectors, index=KO_gene_data.index)


def calculate_module_completeness(
    KeggCalc: keggData.KeggCalculator, KO_gene_data: pd.DataFrame
) -> pd.DataFrame:
    """
    Compute module completeness per bin using NumPy for speed.

    :param KeggCalc: An instance of KeggCalculator containing module definitions.
    :param KO_gene_data: DataFrame containing KO counts per bin with last column "Name".

    :return: DataFrame with completeness values for each module.
    """
    data = KO_gene_data.drop(columns=["Name"]).values
    n_bins = data.shape[0]

    modules = list(KeggCalc.module_definitions.keys())
    n_modules = len(modules)

    # Map KO names to column indices
    col_map = {
        ko: idx for idx, ko in enumerate(KO_gene_data.drop(columns=["Name"]).columns)
    }

    # Prepare result array
    result = np.zeros((n_bins, n_modules), dtype=float)

    for m_idx, module in enumerate(modules):
        # Only keep KOs that exist in the DataFrame

        module_kos = [ko for ko in KeggCalc.module_definitions[module] if ko in col_map]
        if not module_kos:
            continue
        cols = [col_map[ko] for ko in module_kos]

        vals = data[:, cols]
        # vals[vals > 1] = 1  # presence/absence

        result[:, m_idx] = vals.sum(axis=1) / len(KeggCalc.module_definitions[module])

    return pd.DataFrame(result, columns=modules, index=KO_gene_data.index)


def prepare_contig_sizes(contig_to_size: Dict[int, int]) -> np.ndarray:
    """
    Prepare a numpy array of contig sizes for fast access.

    :param contig_to_size: Dictionary mapping contig IDs to contig sizes.

    :return: Numpy array where the index corresponds to the contig ID
             and the value is the contig size.
    """
    max_id = max(contig_to_size)
    contig_sizes = np.zeros(max_id + 1, dtype=np.int64)
    for contig_id, size in contig_to_size.items():
        contig_sizes[contig_id] = size
    return contig_sizes


def compute_N50(lengths: np.ndarray) -> int:
    """
    Compute the N50 value for a given set of contig lengths.

    :param lengths: Numpy array of contig lengths.

    :return: N50 value (contig length at which 50% of the genome is covered).
    """
    arr = np.sort(lengths)
    half = arr.sum() / 2
    csum = np.cumsum(arr)
    return arr[np.searchsorted(csum, half)]


def add_bin_size_and_N50(bins: Iterable[Bin], contig_to_size: Dict[int, int]):
    """
    Add bin size and N50 metrics to a list of bin objects.

    :param bins: List of bin objects.
    :param contig_to_size: Dictionary mapping contig IDs to contig sizes.

    :return: None. The bin objects are updated in place with size and N50.
    """
    # TODO use numpy array everywhere instead of contig_to_size
    contig_sizes = prepare_contig_sizes(contig_to_size)

    for bin_obj in bins:
        lengths = contig_sizes[list(bin_obj.contigs)]  # fast bulk lookup
        total_len = lengths.sum()
        n50 = compute_N50(lengths)

        bin_obj.add_length(int(total_len))
        bin_obj.add_N50(int(n50))


def add_bin_metrics(
    bins: List[Bin],
    contig_info: Dict,
    contamination_weight: float,
    threads: int = 1,
    chunk_size: int = 5000,
):
    """
    Add metrics to a Set of bins.

    :param bins: Set of bin objects.
    :param contig_info: Dictionary containing contig information.
    :param contamination_weight: Weight for contamination assessment.
    :param threads: Number of threads for parallel processing (default is 1).
    :param chunk_size: Number of bins to process in each chunk (default is 5000).

    :return: List of processed bin objects.
    """
    modelPostprocessing = get_modelPostprocessing()
    postProcessor = modelPostprocessing.modelProcessor(threads)

    contig_to_kegg_counter = contig_info["contig_to_kegg_counter"]
    contig_to_cds_count = contig_info["contig_to_cds_count"]
    contig_to_aa_counter = contig_info["contig_to_aa_counter"]
    contig_to_aa_length = contig_info["contig_to_aa_length"]

    logging.info(f"Assessing bin quality for {len(bins)} bins.")
    assess_bins_quality_by_chunk(
        bins,
        contig_to_kegg_counter,
        contig_to_cds_count,
        contig_to_aa_counter,
        contig_to_aa_length,
        contamination_weight,
        postProcessor,
        chunk_size=chunk_size,
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


def assess_bins_quality_by_chunk(
    bins: Iterable[Bin],
    contig_to_kegg_counter: Dict,
    contig_to_cds_count: Dict,
    contig_to_aa_counter: Dict,
    contig_to_aa_length: Dict,
    contamination_weight: float,
    postProcessor=None,
    threads: int = 1,
    chunk_size: int = 2500,
    disable_bar=False,
):
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
    with tqdm(total=len(bins), unit="bin", disable=disable_bar) as pbar:
        for i, chunk_bins_iter in enumerate(chunks(bins, chunk_size)):
            chunk_bins = list(chunk_bins_iter)
            logging.debug(f"chunk {i}: assessing quality of {len(chunk_bins)} bins")
            bins_scored = assess_bins_quality(
                bins=chunk_bins,
                contig_to_kegg_counter=contig_to_kegg_counter,
                contig_to_cds_count=contig_to_cds_count,
                contig_to_aa_counter=contig_to_aa_counter,
                contig_to_aa_length=contig_to_aa_length,
                contamination_weight=contamination_weight,
                postProcessor=postProcessor,
                threads=threads,
            )
            pbar.update(len(bins_scored))


def assess_bins_quality(
    bins: Iterable[Bin],
    contig_to_kegg_counter: Dict,
    contig_to_cds_count: Dict,
    contig_to_aa_counter: Dict,
    contig_to_aa_length: Dict,
    contamination_weight: float,
    postProcessor=None,
    threads: int = 1,
):
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
        modelPostprocessing = get_modelPostprocessing()
        postProcessor = modelPostprocessing.modelProcessor(threads)

    metadata_df = get_bins_metadata_df(
        bins, contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length
    )

    diamond_complete_results, ko_list_length = get_diamond_feature_per_bin_df(
        bins, contig_to_kegg_counter
    )
    diamond_complete_results = diamond_complete_results.drop(columns=["Name"])

    feature_vectors = pd.concat([metadata_df, diamond_complete_results], axis=1)
    feature_vectors = feature_vectors.sort_values(by="Name")

    # 4: Call general model & specific models and derive predictions"""
    modelProcessing = get_modelProcessing()
    modelProc = modelProcessing.modelProcessor(threads)

    vector_array = feature_vectors.iloc[:, 1:].values.astype(float)

    logging.debug("Predicting completeness and contamination using the general model.")
    general_results_comp, general_results_cont = modelProc.run_prediction_general(
        vector_array
    )

    logging.debug("Predicting completeness using the specific model.")
    specific_model_vector_len = (ko_list_length + len(metadata_df.columns)) - 1

    # also retrieve scaled data for CSM calculations
    specific_results_comp, scaled_features = modelProc.run_prediction_specific(
        vector_array, specific_model_vector_len
    )

    logging.debug(
        "Using cosine similarity to reference data to select an appropriate predictor model."
    )

    final_comp, final_cont, models_chosen, csm_array = (
        postProcessor.calculate_general_specific_ratio(
            vector_array[:, 20],
            scaled_features,
            general_results_comp,
            general_results_cont,
            specific_results_comp,
        )
    )

    final_results = feature_vectors[["Name"]].copy()
    final_results["Completeness"] = np.round(final_comp, 2)
    final_results["Contamination"] = np.round(final_cont, 2)

    for bin_obj in bins:
        completeness = final_results.at[bin_obj.contigs_key, "Completeness"]
        contamination = final_results.at[bin_obj.contigs_key, "Contamination"]

        bin_obj.add_quality(completeness, contamination, contamination_weight)

    return bins
