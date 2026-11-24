#!/usr/bin/env python3
import concurrent.futures as cf
import logging
import multiprocessing
import os
from collections import Counter
from collections.abc import Generator, Iterable, Iterator, Sequence
from contextlib import closing, contextmanager
from itertools import chain, islice

import numpy as np
import pandas as pd
from checkm2 import keggData
from pyroaring import BitMap
from rich.progress import Progress

from binette.bin_manager import Bin

logger = logging.getLogger(__name__)

# Suppress unnecessary TensorFlow warnings
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
logging.getLogger("tensorflow").setLevel(logging.FATAL)

# Lazy loaders for checkm2 components that import keras
# These will only be imported when explicitly called
_modelPostprocessing = None
_modelProcessing = None
_keras_initialized = False


def _initialize_keras_environment():
    """Initialize TensorFlow/Keras to ensure thread safety and memory management"""
    global _keras_initialized
    if not _keras_initialized:
        os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"  # Suppress TF warnings

        try:
            # Only import keras-related modules when needed
            import tensorflow as tf

            # Use a single thread for predictions to avoid thread contention
            tf.config.threading.set_intra_op_parallelism_threads(1)
            tf.config.threading.set_inter_op_parallelism_threads(1)

            _keras_initialized = True
            logger.debug("TensorFlow/Keras environment initialized")
        except Exception as e:
            logger.warning(f"Failed to fully initialize TensorFlow environment: {e!s}")


def get_modelPostprocessing():
    """Lazy load modelPostprocessing module only when needed"""
    global _modelPostprocessing
    if _modelPostprocessing is None:
        # Initialize Keras environment
        _initialize_keras_environment()

        # Only import keras when absolutely needed
        from checkm2 import modelPostprocessing

        _modelPostprocessing = modelPostprocessing
    return _modelPostprocessing


def get_modelProcessing():
    """Lazy load modelProcessing module only when needed"""
    global _modelProcessing
    if _modelProcessing is None:
        # Initialize Keras environment
        _initialize_keras_environment()

        # Only import keras when absolutely needed
        from checkm2 import modelProcessing

        _modelProcessing = modelProcessing

    return _modelProcessing


def get_bins_metadata_df(
    contig_sets: Sequence[BitMap],
    contig_to_cds_count: dict[int, int],
    contig_to_aa_counter: dict[int, Counter],
    contig_to_aa_length: dict[int, int],
    kegg_calc: keggData.KeggCalculator | None = None,
) -> pd.DataFrame:
    """
    Generate a DataFrame containing metadata (CDS count, amino acid length and
    composition) for a list of contig sets.

    :param contig_sets: Contig sets (one per bin) to describe.
    :param contig_to_cds_count: Number of CDS per contig.
    :param contig_to_aa_counter: Amino acid composition per contig.
    :param contig_to_aa_length: Total amino acid length per contig.
    :param kegg_calc: KeggCalculator to reuse (one is created when not given).

    :return: DataFrame with one row per contig set, in input order. The "Name" column
        and the index hold the position of the contig set in ``contig_sets``.
    """
    if kegg_calc is None:
        kegg_calc = keggData.KeggCalculator()
    metadata_order = kegg_calc.return_proper_order("Metadata")

    rows = []
    for i, contigs in enumerate(contig_sets):
        cds = 0
        aa_length = 0
        aa_counter = Counter()
        for contig in contigs:
            cds += contig_to_cds_count.get(contig, 0)
            aa_length += contig_to_aa_length.get(contig, 0)
            counter = contig_to_aa_counter.get(contig)
            if counter:
                aa_counter.update(counter)
        row = {"Name": i, "CDS": cds, "AALength": aa_length}
        row.update(aa_counter)
        rows.append(row)

    metadata_df = pd.DataFrame(rows).fillna(0)

    # Ensure column order
    all_cols = ["Name"] + metadata_order
    for col in metadata_order:
        if col not in metadata_df.columns:
            metadata_df[col] = 0

    metadata_df = metadata_df[all_cols].astype(dict.fromkeys(metadata_order, int))
    metadata_df = metadata_df.set_index("Name", drop=False)

    return metadata_df


def get_diamond_feature_per_bin_df(
    contig_sets: Sequence[BitMap],
    contig_to_kegg_counter: dict[int, Counter],
    kegg_calc: keggData.KeggCalculator | None = None,
) -> tuple[pd.DataFrame, int]:
    """
    Generate a DataFrame containing Diamond feature counts per contig set,
    including KEGG KO counts and completeness information for pathways, categories, and modules.

    :param contig_sets: Contig sets (one per bin) to describe.
    :param contig_to_kegg_counter: KO counts per contig.
    :param kegg_calc: KeggCalculator to reuse (one is created when not given).

    :return: DataFrame with one row per contig set, indexed by the position of the
        contig set in ``contig_sets``, and the number of default KOs.
    """
    if kegg_calc is None:
        kegg_calc = keggData.KeggCalculator()
    defaultKOs = kegg_calc.return_default_values_from_category("KO_Genes")
    row_ids = range(len(contig_sets))

    # --- Aggregate KO counters per contig set ---
    row_to_ko_counter = {}
    for i, contigs in enumerate(contig_sets):
        ko_counter = Counter()
        for contig in contigs:
            contig_ko_counter = contig_to_kegg_counter.get(contig)
            if contig_ko_counter:
                ko_counter.update(contig_ko_counter)
        row_to_ko_counter[i] = ko_counter

    # --- Build KO count DataFrame directly ---
    ko_count_per_bin_df = (
        pd.DataFrame.from_dict(row_to_ko_counter, orient="index")
        .reindex(row_ids)  # keep input order
        .fillna(0)
        .astype(int)
    )

    # Ensure all defaultKOs exist
    ko_count_per_bin_df = ko_count_per_bin_df.reindex(
        columns=list(defaultKOs), fill_value=0
    )

    ko_count_per_bin_df["Name"] = ko_count_per_bin_df.index

    # --- Calculate higher-level completeness ---
    logger.debug("Calculating completeness of pathways, categories, and modules")
    KO_pathways = calculate_KO_group(kegg_calc, "KO_Pathways", ko_count_per_bin_df)
    KO_categories = calculate_KO_group(kegg_calc, "KO_Categories", ko_count_per_bin_df)

    KO_modules = calculate_module_completeness(kegg_calc, ko_count_per_bin_df)

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


def prepare_contig_sizes(contig_to_size: dict[int, int]) -> np.ndarray:
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


def add_bin_size_and_N50(bins: Iterable[Bin], contig_to_size: dict[int, int]):
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


def add_bin_coding_density(
    bins: list[Bin], contig_to_coding_length: dict[int, int]
) -> float | None:
    """
    Calculate the coding density of the given bins.

    :param contig_to_coding_length: A dictionary mapping contig IDs to their total coding lengths.

    :return: The coding density of the bin, or None if the length is not set or is zero.
    """
    for bin_obj in bins:
        bin_obj.add_coding_density(contig_to_coding_length)


QualityScores = tuple[np.ndarray, np.ndarray, np.ndarray]
"""Completeness, contamination and CheckM2 model chosen, one value per contig set."""


class QualityScorer:
    """
    Predict completeness and contamination of contig sets with the CheckM2 models.

    The models and the contig annotations are loaded once, when the scorer is created,
    and reused for every contig set scored.
    """

    CONTIG_INFO_KEYS = (
        "contig_to_kegg_counter",
        "contig_to_cds_count",
        "contig_to_aa_counter",
        "contig_to_aa_length",
    )

    def __init__(self, contig_info: dict, threads: int = 1, batch_size: int = 500):
        """
        :param contig_info: Contig annotations: "contig_to_kegg_counter",
            "contig_to_cds_count", "contig_to_aa_counter" and "contig_to_aa_length".
        :param threads: Number of threads used by the CheckM2 models.
        :param batch_size: Maximum number of contig sets sent to the models at once,
            to bound the memory used by feature vectors.
        """
        self.contig_to_kegg_counter = contig_info["contig_to_kegg_counter"]
        self.contig_to_cds_count = contig_info["contig_to_cds_count"]
        self.contig_to_aa_counter = contig_info["contig_to_aa_counter"]
        self.contig_to_aa_length = contig_info["contig_to_aa_length"]
        self.batch_size = batch_size

        self.kegg_calc = keggData.KeggCalculator()
        self.model_processor = get_modelProcessing().modelProcessor(threads)
        self.post_processor = get_modelPostprocessing().modelProcessor(threads)

    def score(self, contig_sets: Sequence[BitMap]) -> QualityScores:
        """
        Score contig sets.

        :param contig_sets: Contig sets (one per bin) to score.
        :return: Completeness, contamination and CheckM2 model chosen, in input order.
        """
        results = [
            self._score_batch(contig_sets[start : start + self.batch_size])
            for start in range(0, len(contig_sets), self.batch_size)
        ]
        if not results:
            empty = np.array([])
            return empty, empty, np.array([], dtype=object)
        completeness, contamination, models = zip(*results, strict=True)
        return (
            np.concatenate(completeness),
            np.concatenate(contamination),
            np.concatenate(models),
        )

    def _score_batch(self, contig_sets: Sequence[BitMap]) -> QualityScores:
        """Score one batch of contig sets. This code is adapted from CheckM2."""
        metadata_df = get_bins_metadata_df(
            contig_sets,
            self.contig_to_cds_count,
            self.contig_to_aa_counter,
            self.contig_to_aa_length,
            self.kegg_calc,
        )

        diamond_complete_results, ko_list_length = get_diamond_feature_per_bin_df(
            contig_sets, self.contig_to_kegg_counter, self.kegg_calc
        )
        diamond_complete_results = diamond_complete_results.drop(columns=["Name"])

        # Both tables have one row per contig set, in input order.
        feature_vectors = pd.concat([metadata_df, diamond_complete_results], axis=1)

        vector_array = feature_vectors.iloc[:, 1:].values.astype(float)

        logger.debug(
            "Predicting completeness and contamination using the general model"
        )
        general_results_comp, general_results_cont = (
            self.model_processor.run_prediction_general(vector_array)
        )

        logger.debug("Predicting completeness using the specific model")
        specific_model_vector_len = (ko_list_length + len(metadata_df.columns)) - 1

        # also retrieve scaled data for CSM calculations
        specific_results_comp, scaled_features = (
            self.model_processor.run_prediction_specific(
                vector_array, specific_model_vector_len
            )
        )

        logger.debug(
            "Using cosine similarity to reference data to select an appropriate predictor model."
        )
        final_comp, final_cont, models_chosen, _ = (
            self.post_processor.calculate_general_specific_ratio(
                vector_array[:, 20],
                scaled_features,
                general_results_comp,
                general_results_cont,
                specific_results_comp,
            )
        )

        return (
            np.round(np.asarray(final_comp, dtype=float), 2),
            np.round(np.asarray(final_cont, dtype=float), 2),
            np.asarray(models_chosen, dtype=object),
        )


# Scorer of a worker process, created once by the pool initializer.
_worker_scorer: QualityScorer | None = None

# Environment variables limiting the threads of numerical libraries in worker processes:
# each worker uses one thread, parallelism comes from the number of workers.
_SINGLE_THREAD_ENV = {
    var: "1"
    for var in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    )
}


def _init_scoring_worker(contig_info: dict, batch_size: int) -> None:
    """Initialize a scoring worker process: load the CheckM2 models once."""
    global _worker_scorer
    _initialize_keras_environment()
    _worker_scorer = QualityScorer(contig_info, threads=1, batch_size=batch_size)


def _score_in_worker(contig_sets: Sequence[BitMap]) -> QualityScores:
    """Score contig sets in a worker process."""
    return _worker_scorer.score(contig_sets)


@contextmanager
def _single_thread_env_for_children() -> Generator[None]:
    """Set single-thread environment variables while worker processes are started."""
    previous = {var: os.environ.get(var) for var in _SINGLE_THREAD_ENV}
    os.environ.update(_SINGLE_THREAD_ENV)
    try:
        yield
    finally:
        for var, value in previous.items():
            if value is None:
                os.environ.pop(var, None)
            else:
                os.environ[var] = value


def score_contig_set_batches(
    batches: Iterable[Sequence[BitMap]],
    contig_info: dict,
    threads: int,
    batch_size: int = 500,
    max_pending: int | None = None,
) -> Iterator[tuple[int, QualityScores]]:
    """
    Score batches of contig sets in worker processes.

    Each worker loads the CheckM2 models and the contig annotations once. Only contig
    sets are sent to workers and only scores are sent back. Batches are read lazily
    from ``batches``: at most ``max_pending`` batches are being scored at any time,
    so a generator can feed batches without all of them being held in memory.

    :param batches: Batches of contig sets to score.
    :param contig_info: Contig annotations (see QualityScorer).
    :param threads: Number of worker processes.
    :param batch_size: Maximum number of contig sets sent to the models at once.
    :param max_pending: Maximum number of batches submitted and not yet returned
        (default: twice the number of workers).

    :return: Iterator of (batch index, scores), in completion order.
    """
    if max_pending is None:
        max_pending = 2 * threads

    batch_iter = enumerate(batches)
    # Workers only need the annotations used by the scorer.
    scorer_info = {key: contig_info[key] for key in QualityScorer.CONTIG_INFO_KEYS}
    # spawn: the parent process may have TensorFlow loaded, which is not fork-safe.
    executor = cf.ProcessPoolExecutor(
        max_workers=threads,
        mp_context=multiprocessing.get_context("spawn"),
        initializer=_init_scoring_worker,
        initargs=(scorer_info, batch_size),
    )
    pending: dict[cf.Future, int] = {}

    def submit_next() -> bool:
        next_batch = next(batch_iter, None)
        if next_batch is None:
            return False
        index, batch = next_batch
        pending[executor.submit(_score_in_worker, batch)] = index
        return True

    try:
        # Worker processes are started by the first submissions (one per submission
        # until all are started): they inherit single-thread settings, the parent
        # process keeps its own.
        with _single_thread_env_for_children():
            while len(pending) < max_pending and submit_next():
                pass

        while pending:
            done, _ = cf.wait(pending, return_when=cf.FIRST_COMPLETED)
            for future in done:
                index = pending.pop(future)
                yield index, future.result()
                submit_next()
    except BaseException:
        # On error (or if the caller stops early), don't score the queued batches.
        executor.shutdown(wait=True, cancel_futures=True)
        raise
    executor.shutdown(wait=True)


def iter_scored_batches(
    contig_sets: Iterable[BitMap],
    contig_info: dict,
    threads: int = 1,
    checkm2_batch_size: int = 500,
    task_size: int | None = None,
    min_bins_for_parallel: int = 6000,
    min_bins_per_worker: int = 3000,
) -> Iterator[tuple[int, list[BitMap], QualityScores]]:
    """
    Score a stream of contig sets in batches.

    Contig sets are read lazily: only the batches being scored are held in memory.
    The first contig sets are read ahead to choose how to score them: in the main
    process when there are few, otherwise in worker processes.

    :param contig_sets: Contig sets to score.
    :param contig_info: Contig annotations (see QualityScorer).
    :param threads: Number of threads. With more than one thread and more than
        ``min_bins_for_parallel`` contig sets, sets are scored in up to ``threads``
        worker processes.
    :param checkm2_batch_size: Maximum number of contig sets sent to the CheckM2 models
        at once, to bound the memory used by feature vectors.
    :param task_size: Number of contig sets per batch. By default, about 4 batches per
        worker, between 500 and 2000 contig sets each, to balance the load.
    :param min_bins_for_parallel: Up to this number of contig sets, sets are scored in
        the main process, avoiding the cost of starting worker processes.
    :param min_bins_per_worker: Each worker process gets at least this many contig
        sets: a worker loads its own copy of the CheckM2 models, which costs memory
        and start-up time, so few sets are scored with fewer workers.

    :return: Iterator of (position of the first contig set of the batch in the input,
        batch of contig sets, scores), in completion order.
    """
    contig_set_iter = iter(contig_sets)
    # Read enough contig sets to know whether there are few of them.
    look_ahead_size = max(min_bins_for_parallel, threads * min_bins_per_worker)
    first_sets = list(islice(contig_set_iter, look_ahead_size + 1))
    # Number of contig sets, when all of them have been read.
    n_sets = len(first_sets) if len(first_sets) <= look_ahead_size else None
    all_sets = chain(first_sets, contig_set_iter)
    del first_sets

    def batches_of(size: int) -> Iterator[tuple[int, list[BitMap]]]:
        start = 0
        while batch := list(islice(all_sets, size)):
            yield start, batch
            start += len(batch)

    if threads == 1 or (n_sets is not None and n_sets <= min_bins_for_parallel):
        logger.debug("Scoring bins in the main process")
        scorer = QualityScorer(
            contig_info, threads=threads, batch_size=checkm2_batch_size
        )
        for start, batch in batches_of(task_size or 2000):
            yield start, batch, scorer.score(batch)
        return

    if n_sets is None:
        n_workers = threads
        task_size = task_size or 2000
    else:
        n_workers = max(1, min(threads, n_sets // min_bins_per_worker))
        if task_size is None:
            task_size = min(2000, max(500, n_sets // (4 * n_workers)))
    logger.debug(
        f"Scoring bins in {n_workers} worker processes, {task_size} bins per task"
    )

    # Batches being scored, by batch index.
    pending_batches: dict[int, tuple[int, list[BitMap]]] = {}

    def batches_to_score() -> Iterator[list[BitMap]]:
        for index, (start, batch) in enumerate(batches_of(task_size)):
            pending_batches[index] = (start, batch)
            yield batch

    # closing: on error, queued batches are cancelled right away.
    with closing(
        score_contig_set_batches(
            batches_to_score(), contig_info, n_workers, batch_size=checkm2_batch_size
        )
    ) as results:
        for index, scores in results:
            start, batch = pending_batches.pop(index)
            yield start, batch, scores


def add_bin_metrics(
    bins: Iterable[Bin],
    contig_info: dict,
    contamination_weight: float,
    threads: int = 1,
    checkm2_batch_size: int = 500,
    disable_progress_bar: bool = False,
    **scoring_options,
) -> list[Bin]:
    """
    Add completeness, contamination, score and CheckM2 model to bins.

    Bins are updated in place.

    :param bins: Bins to assess.
    :param contig_info: Contig annotations (see QualityScorer).
    :param contamination_weight: Weight of contamination in the bin score.
    :param threads: Number of threads (see iter_scored_batches).
    :param checkm2_batch_size: Maximum number of bins sent to the CheckM2 models at once,
        to bound the memory used by feature vectors.
    :param disable_progress_bar: Disable the progress bar if True.
    :param scoring_options: Other options of iter_scored_batches.

    :return: The bins, with quality metrics.
    """
    bins_list = list(bins)
    if not bins_list:
        logger.warning("No bins provided for quality assessment")
        return []

    logger.info(
        f"Assessing bin quality for {len(bins_list)} bins using {threads} threads"
    )

    with (
        Progress(disable=disable_progress_bar) as progress,
        closing(
            iter_scored_batches(
                (bin_obj.contigs for bin_obj in bins_list),
                contig_info,
                threads=threads,
                checkm2_batch_size=checkm2_batch_size,
                **scoring_options,
            )
        ) as results,
    ):
        task = progress.add_task("Assessing bin quality", total=len(bins_list))
        for start, batch, scores in results:
            for bin_obj, completeness, contamination, model in zip(
                bins_list[start : start + len(batch)], *scores, strict=True
            ):
                bin_obj.add_quality(completeness, contamination, contamination_weight)
                bin_obj.add_model(model)
            progress.update(task, advance=len(batch))

    return bins_list


def assess_and_keep_selectable_bins(
    contig_sets: Iterable[BitMap],
    contig_info: dict,
    contamination_weight: float,
    min_completeness: float,
    max_contamination: float,
    threads: int = 1,
    checkm2_batch_size: int = 500,
    disable_progress_bar: bool = False,
    **scoring_options,
) -> list[Bin]:
    """
    Score a stream of contig sets and make bins of those that can be selected.

    Only contig sets with completeness >= ``min_completeness`` and contamination <=
    ``max_contamination`` become bins; the others are dropped once scored, so memory
    depends on the number of bins kept, not on the number of contig sets.

    :param contig_sets: Contig sets to score.
    :param contig_info: Contig annotations (see QualityScorer).
    :param contamination_weight: Weight of contamination in the bin score.
    :param min_completeness: Minimum completeness of the bins kept.
    :param max_contamination: Maximum contamination of the bins kept.
    :param threads: Number of threads (see iter_scored_batches).
    :param checkm2_batch_size: Maximum number of contig sets sent to the CheckM2 models
        at once.
    :param disable_progress_bar: Disable the progress bar if True.
    :param scoring_options: Other options of iter_scored_batches.

    :return: Bins kept, with quality metrics.
    """
    logger.info(f"Assessing intermediate bin quality using {threads} threads")
    kept_bins = []
    n_scored = 0
    with (
        Progress(disable=disable_progress_bar) as progress,
        closing(
            iter_scored_batches(
                contig_sets,
                contig_info,
                threads=threads,
                checkm2_batch_size=checkm2_batch_size,
                **scoring_options,
            )
        ) as results,
    ):
        task = progress.add_task("Assessing intermediate bin quality", total=None)
        for _, batch, scores in results:
            for contigs, completeness, contamination, model in zip(
                batch, *scores, strict=True
            ):
                bin_obj = Bin(contigs, is_original=False)
                bin_obj.add_quality(completeness, contamination, contamination_weight)
                if bin_obj.is_high_quality(min_completeness, max_contamination):
                    bin_obj.add_model(model)
                    kept_bins.append(bin_obj)
            n_scored += len(batch)
            progress.update(task, advance=len(batch))

    logger.info(
        f"{n_scored} intermediate bins assessed, {len(kept_bins)} with completeness >= "
        f"{min_completeness} and contamination <= {max_contamination} kept for bin "
        "selection"
    )
    return kept_bins


def log_selected_bin_info(
    selected_bins: list[Bin],
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
    logger.debug("High quality bins:")
    for sb in selected_bins:
        if sb.is_high_quality(
            min_completeness=hq_min_completeness, max_contamination=hq_max_conta
        ):
            logger.debug(
                f"  {sb} completeness={sb.completeness}, contamination={sb.contamination}"
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
    logger.info(
        f"{hq_bins}/{len(selected_bins)} selected bins have high quality {thresholds}"
    )
