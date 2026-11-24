#!/usr/bin/env python3
import concurrent.futures as cf
import logging
import multiprocessing
import os
from collections.abc import Generator, Iterable, Iterator, Sequence
from contextlib import closing, contextmanager
from dataclasses import dataclass
from itertools import chain, islice

import numpy as np
import scipy.sparse
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


class FeatureBuilder:
    """
    Build the CheckM2 feature vectors of contig sets.

    The features of a contig set are sums over its contigs (CDS count, amino acid
    length and composition, KO counts), plus the completeness of KEGG pathways,
    modules and categories computed from the KO counts. Per-contig values are stored
    in matrices once, so the features of a batch of contig sets come from a few
    integer matrix products with the contig set x contig incidence matrix.
    """

    def __init__(self, contig_info: dict, kegg_calc: keggData.KeggCalculator):
        """
        :param contig_info: Contig annotations: "contig_to_kegg_counter",
            "contig_to_cds_count", "contig_to_aa_counter" and "contig_to_aa_length".
        :param kegg_calc: CheckM2 KeggCalculator, giving feature names and KEGG groups.
        """
        contig_to_kegg_counter = contig_info["contig_to_kegg_counter"]
        contig_to_cds_count = contig_info["contig_to_cds_count"]
        contig_to_aa_counter = contig_info["contig_to_aa_counter"]
        contig_to_aa_length = contig_info["contig_to_aa_length"]

        self.metadata_names = kegg_calc.return_proper_order("Metadata")
        self.ko_names = list(kegg_calc.return_default_values_from_category("KO_Genes"))
        ko_index = {ko: i for i, ko in enumerate(self.ko_names)}

        # Contigs without annotation are not stored: they add nothing to features.
        self.n_contigs = 1 + max(
            max(annotations, default=-1)
            for annotations in (
                contig_to_kegg_counter,
                contig_to_cds_count,
                contig_to_aa_counter,
                contig_to_aa_length,
            )
        )

        # Contig x metadata (CDS count, amino acid length and composition).
        metadata_index = {name: i for i, name in enumerate(self.metadata_names)}
        aa_names = set(self.metadata_names) - {"CDS", "AALength"}
        self.contig_metadata = np.zeros(
            (self.n_contigs, len(self.metadata_names)), dtype=np.int64
        )
        for contig, count in contig_to_cds_count.items():
            self.contig_metadata[contig, metadata_index["CDS"]] = count
        for contig, length in contig_to_aa_length.items():
            self.contig_metadata[contig, metadata_index["AALength"]] = length
        for contig, aa_counter in contig_to_aa_counter.items():
            for aa, count in aa_counter.items():
                if aa in aa_names:
                    self.contig_metadata[contig, metadata_index[aa]] += count

        # Contig x KO counts (KOs not used by CheckM2 are ignored).
        rows, cols, counts = [], [], []
        for contig, ko_counter in contig_to_kegg_counter.items():
            for ko, count in ko_counter.items():
                if ko in ko_index:
                    rows.append(contig)
                    cols.append(ko_index[ko])
                    counts.append(count)
        self.contig_ko_counts = scipy.sparse.csr_matrix(
            (np.array(counts, dtype=np.int64), (rows, cols)),
            shape=(self.n_contigs, len(self.ko_names)),
        )

        # KO x group membership and group sizes. Pathway and category completeness
        # count the KOs present; module completeness sums KO counts, and a KO listed
        # twice in a module definition counts twice.
        self.pathway_names, self.pathway_members, self.pathway_sizes = (
            self._group_matrix(ko_index, self._group_kos(kegg_calc, "KO_Pathways"))
        )
        self.category_names, self.category_members, self.category_sizes = (
            self._group_matrix(ko_index, self._group_kos(kegg_calc, "KO_Categories"))
        )
        self.module_names, self.module_members, self.module_sizes = self._group_matrix(
            ko_index, kegg_calc.module_definitions
        )

    @staticmethod
    def _group_kos(kegg_calc: keggData.KeggCalculator, group: str) -> dict:
        """KOs of each entry of a KEGG group (pathways or categories), in CheckM2 order."""
        mapping = kegg_calc.path_category_mapping
        kos_by_entry = mapping.groupby(group, sort=False)["Kegg_ID"].agg(list)
        return {
            entry: kos_by_entry.get(entry, [])
            for entry in kegg_calc.return_default_values_from_category(group)
        }

    @staticmethod
    def _group_matrix(
        ko_index: dict[str, int], group_kos: dict[str, list[str]]
    ) -> tuple[list[str], scipy.sparse.csr_matrix, np.ndarray]:
        """
        KO x group matrix (number of times a KO is listed in a group) and group sizes.

        :param ko_index: Column of each KO in the KO count matrix.
        :param group_kos: KOs of each group, in group order.
        :return: Group names, KO x group matrix and number of KOs listed per group.
        """
        names = list(group_kos)
        rows, cols = [], []
        sizes = np.zeros(len(names), dtype=np.int64)
        for j, name in enumerate(names):
            kos = list(group_kos[name])
            sizes[j] = len(kos)
            for ko in kos:
                if ko in ko_index:
                    rows.append(ko_index[ko])
                    cols.append(j)
        members = scipy.sparse.csr_matrix(
            (np.ones(len(rows), dtype=np.int64), (rows, cols)),
            shape=(len(ko_index), len(names)),
        )
        # Empty groups have no members: any size other than 0 gives 0.
        sizes[sizes == 0] = 1
        return names, members, sizes

    @property
    def feature_names(self) -> list[str]:
        """Names of the features, in the order of CheckM2 feature vectors."""
        return (
            self.metadata_names
            + self.ko_names
            + self.pathway_names
            + self.module_names
            + self.category_names
        )

    @property
    def specific_model_vector_len(self) -> int:
        """Number of leading features used by the CheckM2 specific model."""
        return len(self.metadata_names) + len(self.ko_names)

    def incidence_matrix(
        self, contig_sets: Sequence[BitMap]
    ) -> scipy.sparse.csr_matrix:
        """Contig set x contig matrix: 1 where the contig is in the set."""
        contig_ids = [
            np.asarray(contigs.to_array(), dtype=np.int64) for contigs in contig_sets
        ]
        indptr = np.zeros(len(contig_sets) + 1, dtype=np.int64)
        indptr[1:] = np.cumsum([len(ids) for ids in contig_ids])
        indices = (
            np.concatenate(contig_ids) if contig_ids else np.zeros(0, dtype=np.int64)
        )
        # Contigs without annotation are beyond the stored matrices: drop them.
        if indices.size and indices.max() >= self.n_contigs:
            annotated = indices < self.n_contigs
            kept_before = np.concatenate(([0], np.cumsum(annotated)))
            indptr = kept_before[indptr]
            indices = indices[annotated]
        return scipy.sparse.csr_matrix(
            (np.ones(len(indices), dtype=np.int64), indices, indptr),
            shape=(len(contig_sets), self.n_contigs),
        )

    def build(self, contig_sets: Sequence[BitMap]) -> np.ndarray:
        """
        Feature vectors of contig sets.

        :param contig_sets: Contig sets (one per bin).
        :return: Array with one row per contig set, in input order, and one column
            per feature (see feature_names).
        """
        incidence = self.incidence_matrix(contig_sets)
        metadata = incidence @ self.contig_metadata
        ko_counts = incidence @ self.contig_ko_counts
        ko_presence = ko_counts.copy()
        ko_presence.eliminate_zeros()
        ko_presence.data[:] = 1

        pathways = (ko_presence @ self.pathway_members).toarray() / self.pathway_sizes
        modules = (ko_counts @ self.module_members).toarray() / self.module_sizes
        categories = (
            ko_presence @ self.category_members
        ).toarray() / self.category_sizes

        return np.hstack(
            [
                np.asarray(metadata, dtype=float),
                ko_counts.toarray().astype(float),
                pathways,
                modules,
                categories,
            ]
        )


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


@dataclass(frozen=True)
class SelectionThresholds:
    """Completeness and contamination a bin needs to be selectable."""

    min_completeness: float
    max_contamination: float


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

    def __init__(
        self,
        contig_info: dict,
        threads: int = 1,
        batch_size: int = 500,
        selectable_only: SelectionThresholds | None = None,
    ):
        """
        :param contig_info: Contig annotations: "contig_to_kegg_counter",
            "contig_to_cds_count", "contig_to_aa_counter" and "contig_to_aa_length".
        :param threads: Number of threads used by the CheckM2 models.
        :param batch_size: Maximum number of contig sets sent to the models at once,
            to bound the memory used by feature vectors.
        :param selectable_only: When given, the CheckM2 steps that can't change
            whether a contig set is selectable are skipped (see _score_batch). Contig
            sets found not selectable get a NaN completeness and no model.
        """
        self.selectable_only = selectable_only
        self.batch_size = batch_size

        self.features = FeatureBuilder(contig_info, keggData.KeggCalculator())
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
        vector_array = self.features.build(contig_sets)

        logger.debug(
            "Predicting completeness and contamination using the general model"
        )
        general_comp, general_cont = self.model_processor.run_prediction_general(
            vector_array
        )

        n_sets = len(contig_sets)
        completeness = np.full(n_sets, np.nan)
        contamination = np.round(np.asarray(general_cont, dtype=float), 2)
        models = np.full(n_sets, None, dtype=object)

        # Rows still to score. With selectable_only, two shortcuts skip the costly
        # steps for contig sets that can't be selected, without changing any result:
        # - contamination only comes from the general model: above the maximum
        #   contamination, the set is not selectable whatever its completeness;
        # - the final completeness is the general or the specific prediction: when
        #   both are below the minimum completeness, the model choice can't help.
        rows = np.arange(n_sets)
        thresholds = self.selectable_only
        if thresholds is not None:
            rows = rows[contamination <= thresholds.max_contamination]
        if rows.size == 0:
            return completeness, contamination, models

        logger.debug("Predicting completeness using the specific model")
        specific_model_vector_len = self.features.specific_model_vector_len

        # also retrieve scaled data for CSM calculations
        specific_comp, scaled_features = self.model_processor.run_prediction_specific(
            vector_array[rows], specific_model_vector_len
        )

        if thresholds is not None:
            best_comp = np.maximum(general_comp[rows], specific_comp)
            can_pass = (
                np.round(best_comp.astype(float), 2) >= thresholds.min_completeness
            )
            rows = rows[can_pass]
            specific_comp = specific_comp[can_pass]
            scaled_features = scaled_features[can_pass]
            if rows.size == 0:
                return completeness, contamination, models

        logger.debug(
            "Using cosine similarity to reference data to select an appropriate predictor model."
        )
        final_comp, _, models_chosen, _ = (
            self.post_processor.calculate_general_specific_ratio(
                vector_array[rows, 20],
                scaled_features,
                general_comp[rows],
                general_cont[rows],
                specific_comp,
            )
        )

        completeness[rows] = np.round(np.asarray(final_comp, dtype=float), 2)
        models[rows] = models_chosen
        return completeness, contamination, models


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


def _init_scoring_worker(
    contig_info: dict, batch_size: int, selectable_only: SelectionThresholds | None
) -> None:
    """Initialize a scoring worker process: load the CheckM2 models once."""
    global _worker_scorer
    _initialize_keras_environment()
    _worker_scorer = QualityScorer(
        contig_info,
        threads=1,
        batch_size=batch_size,
        selectable_only=selectable_only,
    )


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
    selectable_only: SelectionThresholds | None = None,
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
    :param selectable_only: See QualityScorer.

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
        initargs=(scorer_info, batch_size, selectable_only),
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
    min_bins_per_worker: int = 1000,
    selectable_only: SelectionThresholds | None = None,
) -> Iterator[tuple[int, list[BitMap], QualityScores]]:
    """
    Score a stream of contig sets in batches.

    Contig sets are read lazily: only the batches being scored are held in memory.
    The first contig sets are read ahead to choose how to score them: in the main
    process when there are few, otherwise in worker processes.

    :param contig_sets: Contig sets to score.
    :param contig_info: Contig annotations (see QualityScorer).
    :param threads: Maximum number of worker processes.
    :param checkm2_batch_size: Maximum number of contig sets sent to the CheckM2 models
        at once, to bound the memory used by feature vectors.
    :param task_size: Number of contig sets per batch. By default, about 4 batches per
        worker, between 500 and 2000 contig sets each, to balance the load.
    :param min_bins_per_worker: Each worker process gets at least this many contig
        sets: a worker takes a few seconds to start and loads its own copy of the
        CheckM2 models (about 1 GB), so few sets are scored with fewer workers. When
        that leaves fewer than two workers, sets are scored in the main process.
    :param selectable_only: See QualityScorer.

    :return: Iterator of (position of the first contig set of the batch in the input,
        batch of contig sets, scores), in completion order.
    """
    contig_set_iter = iter(contig_sets)
    # Read enough contig sets to know whether there are too few for all workers.
    look_ahead_size = threads * min_bins_per_worker
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

    if n_sets is None:
        n_workers = threads
    else:
        n_workers = min(threads, n_sets // min_bins_per_worker)

    if n_workers < 2:
        logger.debug("Scoring bins in the main process")
        scorer = QualityScorer(
            contig_info,
            threads=threads,
            batch_size=checkm2_batch_size,
            selectable_only=selectable_only,
        )
        for start, batch in batches_of(task_size or 2000):
            yield start, batch, scorer.score(batch)
        return

    if task_size is None:
        task_size = (
            2000 if n_sets is None else min(2000, max(500, n_sets // (4 * n_workers)))
        )
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
            batches_to_score(),
            contig_info,
            n_workers,
            batch_size=checkm2_batch_size,
            selectable_only=selectable_only,
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
    # Contig sets found not selectable before the end of the CheckM2 scoring.
    n_set_aside_early = 0
    with (
        Progress(disable=disable_progress_bar) as progress,
        closing(
            iter_scored_batches(
                contig_sets,
                contig_info,
                threads=threads,
                checkm2_batch_size=checkm2_batch_size,
                selectable_only=SelectionThresholds(
                    min_completeness, max_contamination
                ),
                **scoring_options,
            )
        ) as results,
    ):
        task = progress.add_task("Assessing intermediate bin quality", total=None)
        for _, batch, scores in results:
            for contigs, completeness, contamination, model in zip(
                batch, *scores, strict=True
            ):
                if np.isnan(completeness):
                    n_set_aside_early += 1
                    continue
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
    logger.debug(
        f"{n_set_aside_early} intermediate bins found not selectable before the "
        "end of the CheckM2 scoring"
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
