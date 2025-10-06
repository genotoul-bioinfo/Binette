from collections import Counter

from pyroaring import BitMap

from binette import bin_quality
from binette.bin_manager import Bin
from binette.bin_quality import (
    balanced_chunks,
    chunks,
)


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


def test_balanced_chunks_normal_distribution():
    """Test balanced_chunks with normal distribution scenarios."""

    # Test case: 50 bins with 8 threads
    bins_50 = list(range(50))
    chunks_50_8 = list(balanced_chunks(bins_50, 8))

    assert len(chunks_50_8) == 8
    # 50 / 8 = 6 remainder 2, so first 2 chunks get 7 items, rest get 6
    chunk_sizes = [len(chunk) for chunk in chunks_50_8]
    assert chunk_sizes == [7, 7, 6, 6, 6, 6, 6, 6]
    # Verify all items are included
    all_items = [item for chunk in chunks_50_8 for item in chunk]
    assert sorted(all_items) == bins_50


def test_balanced_chunks_large_dataset():
    """Test balanced_chunks with larger datasets."""

    # Test case: 100 bins with 8 threads
    bins_100 = list(range(100))
    chunks_100_8 = list(balanced_chunks(bins_100, 8))

    assert len(chunks_100_8) == 8
    # 100 / 8 = 12 remainder 4, so first 4 chunks get 13 items, rest get 12
    chunk_sizes = [len(chunk) for chunk in chunks_100_8]
    assert chunk_sizes == [13, 13, 13, 13, 12, 12, 12, 12]
    # Verify all items are included
    all_items = [item for chunk in chunks_100_8 for item in chunk]
    assert sorted(all_items) == bins_100


def test_balanced_chunks_fewer_items_than_chunks():
    """Test balanced_chunks when there are fewer items than requested chunks."""

    bins_5 = list(range(5))
    chunks_5_8 = list(balanced_chunks(bins_5, 8))

    assert len(chunks_5_8) == 5  # Should only create 5 chunks
    chunk_sizes = [len(chunk) for chunk in chunks_5_8]
    assert chunk_sizes == [1, 1, 1, 1, 1]
    # Verify all items are included
    all_items = [item for chunk in chunks_5_8 for item in chunk]
    assert sorted(all_items) == bins_5


def test_balanced_chunks_empty_list():
    """Test balanced_chunks with empty input."""

    bins_empty = []
    chunks_empty = list(balanced_chunks(bins_empty, 8))

    assert chunks_empty == []


def test_balanced_chunks_single_item():
    """Test balanced_chunks with single item."""

    bins_1 = [42]
    chunks_1 = list(balanced_chunks(bins_1, 8))

    assert len(chunks_1) == 1
    assert chunks_1[0] == [42]


def test_balanced_chunks_exact_division():
    """Test balanced_chunks when items divide evenly into chunks."""

    bins_16 = list(range(16))
    chunks_16_4 = list(balanced_chunks(bins_16, 4))

    assert len(chunks_16_4) == 4
    chunk_sizes = [len(chunk) for chunk in chunks_16_4]
    assert chunk_sizes == [4, 4, 4, 4]
    # Verify all items are included
    all_items = [item for chunk in chunks_16_4 for item in chunk]
    assert sorted(all_items) == bins_16


def test_balanced_chunks_string_items():
    """Test balanced_chunks with non-numeric items."""

    bins_strings = ["bin_a", "bin_b", "bin_c", "bin_d", "bin_e"]
    chunks_strings = list(balanced_chunks(bins_strings, 3))

    assert len(chunks_strings) == 3
    # 5 / 3 = 1 remainder 2, so first 2 chunks get 2 items, last gets 1
    chunk_sizes = [len(chunk) for chunk in chunks_strings]
    assert chunk_sizes == [2, 2, 1]
    # Verify all items are included
    all_items = [item for chunk in chunks_strings for item in chunk]
    assert sorted(all_items) == sorted(bins_strings)


def test_balanced_chunks_single_thread():
    """Test balanced_chunks with single thread (one chunk)."""

    bins_3 = [1, 2, 3]
    chunks_3_1 = list(balanced_chunks(bins_3, 1))

    assert len(chunks_3_1) == 1
    assert chunks_3_1[0] == [1, 2, 3]


def test_balanced_chunks_distribution_properties():
    """Test that balanced_chunks maintains balanced distribution properties."""

    # Test the specific example from the conversation - 100 bins, 8 threads
    # This should create 8 chunks instead of 9 to match the number of threads
    bins_100_example = list(range(100))
    chunks_100_example = list(balanced_chunks(bins_100_example, 8))

    assert len(chunks_100_example) == 8
    # Verify the distribution is balanced
    chunk_sizes_example = [len(chunk) for chunk in chunks_100_example]
    min_size = min(chunk_sizes_example)
    max_size = max(chunk_sizes_example)
    assert max_size - min_size <= 1  # Difference should be at most 1
    assert sum(chunk_sizes_example) == 100  # All items included
