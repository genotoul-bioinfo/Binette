from calendar import c
from itertools import islice
from binette import bin_quality

from collections import Counter
import pandas as pd
from unittest.mock import Mock, patch

from unittest.mock import Mock, patch
from binette.bin_quality import (
    Bin,
    add_bin_metrics,
    assess_bins_quality_by_chunk,
    assess_bins_quality,
    chunks,
    get_diamond_feature_per_bin_df,
    get_bins_metadata_df,
)
from pyroaring import BitMap
from checkm2 import keggData, modelPostprocessing, modelProcessing

from unittest.mock import Mock, patch, MagicMock


def test_compute_N50():
    assert bin_quality.compute_N50([50]) == 50
    assert bin_quality.compute_N50([0]) == 0
    assert bin_quality.compute_N50([30, 40, 30]) == 30
    assert bin_quality.compute_N50([1, 3, 3, 4, 5, 5, 6, 9, 10, 24]) == 9


def test_chunks():
    # Test case 1
    iterable_1 = [1, 2, 3, 4, 5, 6]
    size_1 = 2
    expected_output_1 = [(1, 2), (3, 4), (5, 6)]

    result_1 = list(chunks(iterable_1, size_1))
    assert result_1 == expected_output_1

    # Test case 2
    iterable_2 = [10, 20, 30, 40, 50]
    size_2 = 3
    expected_output_2 = [(10, 20, 30), (40, 50)]

    result_2 = list(chunks(iterable_2, size_2))
    assert result_2 == expected_output_2

    # Test case 3 (Empty iterable)
    iterable_3 = []
    size_3 = 5
    expected_output_3 = []

    result_3 = list(chunks(iterable_3, size_3))
    assert result_3 == expected_output_3

    # Test case 4 (Iterable length less than chunk size)
    iterable_4 = [100, 200, 300]
    size_4 = 5
    expected_output_4 = [(100, 200, 300)]

    result_4 = list(chunks(iterable_4, size_4))
    assert result_4 == expected_output_4


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

    def add_N50(self, N50):
        self.N50 = N50

    def add_quality(self, comp, cont, weight):

        self.completeness = comp
        self.contamination = cont
        self.score = comp - weight * cont


def test_get_bins_metadata_df():
    # Mock input data
    bins = [Bin(BitMap((1, 3))), Bin(BitMap((2,)))]

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
        bins, contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length
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
    bins = [Bin(BitMap((1, 2))), Bin(BitMap((2, 3)))]

    contig_to_kegg_counter = {
        1: Counter({"K01810": 5, "K15916": 7}),
        2: Counter({"K01810": 10}),
        3: Counter({"K00918": 8}),
    }

    # Call the function
    result_df, default_ko_count = bin_quality.get_diamond_feature_per_bin_df(
        bins, contig_to_kegg_counter
    )

    assert (
        result_df.loc[bins[0].contigs_key, "K01810"] == 15
    )  # in bin1 from contig 1 and 2
    assert result_df.loc[bins[0].contigs_key, "K15916"] == 7  # in bin1 from contig 1
    assert (
        result_df.loc[bins[1].contigs_key, "K01810"] == 10
    )  # this ko is not in any contig of bin 2
    assert result_df.loc[bins[1].contigs_key, "K00918"] == 8  # in bin2 from contig 3


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
