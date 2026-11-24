import logging
from collections import Counter
from unittest.mock import patch

import numpy as np
import pytest
from pyroaring import BitMap

from binette import bin_quality
from binette.bin_manager import Bin


def test_compute_N50():
    assert bin_quality.compute_N50([50]) == 50
    assert bin_quality.compute_N50([0]) == 0
    assert bin_quality.compute_N50([30, 40, 30]) == 30
    assert bin_quality.compute_N50([1, 3, 3, 4, 5, 5, 6, 9, 10, 24]) == 9


class BinOLD:
    def __init__(self, bin_id, contigs):
        self.id = bin_id
        self.contigs = contigs
        self.length = 0  # Mocking the add_length method
        self.N50 = 0  # Mocking the add_N50 method

    def add_length(self, length):
        self.length = length

    def add_N50(self, N50):
        self.N50 = N50

    def add_quality(self, comp, cont, weight):
        self.completeness = comp
        self.contamination = cont
        self.score = comp - weight * cont


def test_get_bins_metadata_df():
    # Mock input data
    contig_sets = [BitMap((1, 3)), BitMap((2,))]

    contig_to_cds_count = {1: 10, 2: 45, 3: 20, 4: 25}
    contig_to_aa_counter = {
        1: Counter({"A": 5, "D": 10}),
        2: Counter({"G": 8, "V": 12, "T": 2}),
        3: Counter({"D": 8, "Y": 12}),
    }
    contig_to_aa_length = {
        1: 1000,
        2: 1500,
        3: 2000,
        4: 2500,
    }

    # Call the function
    result_df = bin_quality.get_bins_metadata_df(
        contig_sets, contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length
    )

    # Define expected values based on the provided input
    expected_columns = [
        "Name",
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
        "AALength",
        "CDS",
    ]

    expected_values = [
        ["NA", 5, 0, 18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 3000, 30],
        ["NA", 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 12, 0, 0, 1500, 45],
    ]

    result_df["Name"] = "NA"
    result_df.index = range(len(result_df))
    print(result_df)
    # Check if the generated DataFrame matches the expected DataFrame
    assert result_df.columns.tolist() == expected_columns
    assert result_df.values.tolist() == expected_values


def test_get_diamond_feature_per_bin_df():
    # Mock input data
    contig_sets = [BitMap((1, 2)), BitMap((2, 3))]

    contig_to_kegg_counter = {
        1: Counter({"K01810": 5, "K15916": 7}),
        2: Counter({"K01810": 10}),
        3: Counter({"K00918": 8}),
    }

    # Call the function
    result_df, default_ko_count = bin_quality.get_diamond_feature_per_bin_df(
        contig_sets, contig_to_kegg_counter
    )

    # rows are indexed by the position of the contig set
    assert result_df.index.tolist() == [0, 1]
    assert result_df.loc[0, "K01810"] == 15  # in set 0 from contig 1 and 2
    assert result_df.loc[0, "K15916"] == 7  # in set 0 from contig 1
    assert result_df.loc[1, "K01810"] == 10  # in set 1 from contig 2
    assert result_df.loc[1, "K00918"] == 8  # in set 1 from contig 3


def test_add_bin_size_and_N50():
    # Mock input data
    bins = [Bin(BitMap((1, 2))), Bin(BitMap((2, 3)))]

    contig_to_size = {
        1: 1000,
        2: 1500,
        3: 2000,
    }

    # Call the function
    bin_quality.add_bin_size_and_N50(bins, contig_to_size)

    # Assertions to verify if add_length and add_N50 were called with the correct values
    assert bins[0].length == 2500
    assert bins[0].N50 == 1500
    assert bins[1].length == 3500
    assert bins[1].N50 == 2000


class FakeScorer:
    """Stand-in for QualityScorer: scores derived from the contig ids, no CheckM2."""

    def __init__(self, contig_info, threads=1, batch_size=500, selectable_only=None):
        pass

    def score(self, contig_sets):
        first = np.array([min(c) for c in contig_sets], dtype=float)
        models = np.array([Bin.CHECKM2_MODELS[int(f) % 2] for f in first], dtype=object)
        return first, first / 10, models


def test_add_bin_metrics_sequential_path_assigns_each_bin_its_scores():
    bins = [Bin(BitMap((i,))) for i in range(25)]

    with patch("binette.bin_quality.QualityScorer", FakeScorer):
        result_bins = bin_quality.add_bin_metrics(
            bins=bins,
            contig_info={},
            contamination_weight=2,
            threads=1,
            disable_progress_bar=True,
            task_size=4,  # several tasks, the last one partial
        )

    assert result_bins == bins
    for i, bin_obj in enumerate(result_bins):
        assert bin_obj.completeness == i
        assert bin_obj.contamination == i / 10
        assert bin_obj.score == i - 2 * (i / 10)
        assert bin_obj.checkm2_model == Bin.CHECKM2_MODELS[i % 2]


def make_contig_info(n_contigs):
    """Contig annotations with real KO identifiers, so CheckM2 models can run."""
    from checkm2 import keggData

    kos = list(
        keggData.KeggCalculator().return_default_values_from_category("KO_Genes")
    )
    aas = "ACDEFGHIKLMNPQRSTVWY"
    return {
        "contig_to_kegg_counter": {
            c: Counter({kos[(c * 7 + k) % len(kos)]: 1 + k for k in range(c % 4)})
            for c in range(n_contigs)
        },
        "contig_to_cds_count": {c: 1 + c % 5 for c in range(n_contigs)},
        "contig_to_aa_counter": {
            c: Counter({aa: 10 + (c * j) % 50 for j, aa in enumerate(aas)})
            for c in range(n_contigs)
        },
        "contig_to_aa_length": {
            c: sum(10 + (c * j) % 50 for j in range(len(aas))) for c in range(n_contigs)
        },
    }


def make_contig_sets(n_sets, n_contigs):
    return [
        BitMap(range(i % n_contigs, min(n_contigs, i % n_contigs + 5 + i % 40)))
        for i in range(n_sets)
    ]


def test_quality_scorer_results_do_not_depend_on_batches():
    """Real CheckM2 models: batch size and input order do not change scores."""
    contig_info = make_contig_info(300)
    contig_sets = make_contig_sets(60, 300)

    one_batch = bin_quality.QualityScorer(contig_info, batch_size=500).score(
        contig_sets
    )
    small_batches = bin_quality.QualityScorer(contig_info, batch_size=7).score(
        contig_sets
    )
    reversed_input = bin_quality.QualityScorer(contig_info, batch_size=500).score(
        contig_sets[::-1]
    )

    for a, b, r in zip(one_batch, small_batches, reversed_input, strict=True):
        assert len(a) == len(contig_sets)
        assert list(a) == list(b)
        assert list(a) == list(r[::-1])


def test_selectable_only_scoring_keeps_scores_of_selectable_bins():
    """Real CheckM2 models: skipping steps only drops bins that can't be selected."""
    contig_info = make_contig_info(300)
    contig_sets = make_contig_sets(200, 300)
    full_comp, full_cont, full_models = bin_quality.QualityScorer(contig_info).score(
        contig_sets
    )
    # Thresholds in the middle of the scores, so that both shortcuts are used.
    thresholds = bin_quality.SelectionThresholds(
        min_completeness=float(np.median(full_comp)),
        max_contamination=float(np.median(full_cont)),
    )

    comp, cont, models = bin_quality.QualityScorer(
        contig_info, selectable_only=thresholds
    ).score(contig_sets)

    selectable = (full_comp >= thresholds.min_completeness) & (
        full_cont <= thresholds.max_contamination
    )
    set_aside = np.isnan(comp)
    assert selectable.any() and set_aside.any()
    assert not (selectable & set_aside).any()
    too_contaminated = full_cont > thresholds.max_contamination
    assert (set_aside & too_contaminated).any()  # first shortcut used
    assert (set_aside & ~too_contaminated).any()  # second shortcut used
    assert list(cont) == list(full_cont)
    assert list(comp[~set_aside]) == list(full_comp[~set_aside])
    assert list(models[~set_aside]) == list(full_models[~set_aside])
    assert all(m is None for m in models[set_aside])


def test_add_bin_metrics_parallel_path_matches_sequential_path():
    """Real CheckM2 models in worker processes give the same scores as in-process."""
    contig_info = make_contig_info(300)
    contig_sets = make_contig_sets(90, 300)

    sequential = bin_quality.add_bin_metrics(
        [Bin(c) for c in contig_sets],
        contig_info,
        contamination_weight=2,
        threads=1,
        disable_progress_bar=True,
    )
    parallel = bin_quality.add_bin_metrics(
        [Bin(c) for c in contig_sets],
        contig_info,
        contamination_weight=2,
        threads=2,
        disable_progress_bar=True,
        task_size=20,  # 5 tasks, the last one partial
        min_bins_per_worker=1,  # use both workers
    )

    for s, p in zip(sequential, parallel, strict=True):
        assert s.contigs == p.contigs
        assert (s.completeness, s.contamination, s.score, s.checkm2_model) == (
            p.completeness,
            p.contamination,
            p.score,
            p.checkm2_model,
        )


def test_assess_and_keep_selectable_bins_keeps_only_selectable_bins():
    contig_sets = (BitMap((i,)) for i in range(25))  # a generator, read lazily

    with patch("binette.bin_quality.QualityScorer", FakeScorer):
        kept = bin_quality.assess_and_keep_selectable_bins(
            contig_sets,
            contig_info={},
            contamination_weight=2,
            min_completeness=10,
            max_contamination=2,
            threads=1,
            disable_progress_bar=True,
            task_size=4,
        )

    # completeness = i, contamination = i / 10
    assert [b.contigs for b in kept] == [BitMap((i,)) for i in range(10, 21)]
    for bin_obj in kept:
        i = bin_obj.contigs.min()
        assert (bin_obj.completeness, bin_obj.contamination) == (i, i / 10)
        assert bin_obj.checkm2_model == Bin.CHECKM2_MODELS[i % 2]
        assert not bin_obj.is_original


def test_assess_and_keep_selectable_bins_streams_to_worker_processes():
    """More contig sets than read ahead: scored by workers as they are read."""
    contig_info = make_contig_info(300)
    contig_sets = make_contig_sets(90, 300)

    sequential = bin_quality.add_bin_metrics(
        [Bin(c) for c in contig_sets],
        contig_info,
        contamination_weight=2,
        threads=1,
        disable_progress_bar=True,
    )
    streamed = bin_quality.assess_and_keep_selectable_bins(
        iter(contig_sets),
        contig_info,
        contamination_weight=2,
        min_completeness=0,
        max_contamination=100,  # keep all bins
        threads=2,
        disable_progress_bar=True,
        task_size=20,
        min_bins_per_worker=10,  # read ahead 20 contig sets only
    )

    # Bins come back in batch completion order.
    streamed_by_contigs = {b.contigs_key: b for b in streamed}
    assert len(streamed_by_contigs) == len(sequential)
    for s in sequential:
        p = streamed_by_contigs[s.contigs_key]
        assert (s.completeness, s.contamination, s.score, s.checkm2_model) == (
            p.completeness,
            p.contamination,
            p.score,
            p.checkm2_model,
        )


def test_add_bin_metrics_empty_bins():
    result_bins = bin_quality.add_bin_metrics(
        bins=[],
        contig_info=[],
        contamination_weight=2,
        threads=1,  # Single thread
        checkm2_batch_size=500,
        disable_progress_bar=True,
    )

    assert result_bins == [], "Result should be an empty list when input bins are empty"


@pytest.fixture
def bins():
    b1 = Bin(contigs=BitMap({1}), origin="set1", name="bin1")
    b2 = Bin(contigs=BitMap({3}), origin="set1", name="bin2")
    b3 = Bin(contigs=BitMap({3, 2}), origin="set1", name="bin3")

    b1.add_quality(100, 0, 0)
    b2.add_quality(95, 10, 0)
    b3.add_quality(70, 20, 0)

    return [b1, b2, b3]


def test_log_selected_bin_info(caplog, bins):
    caplog.set_level(logging.INFO)

    hq_min_completeness = 85
    hq_max_conta = 15

    # Call the function
    bin_quality.log_selected_bin_info(bins, hq_min_completeness, hq_max_conta)

    # Check if the logs contain expected messages
    expected_logs = "2/3 selected bins have high quality (completeness >= 85 and contamination <= 15)"

    assert expected_logs in caplog.text
