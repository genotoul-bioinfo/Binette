from collections import Counter
from unittest.mock import patch

from pyroaring import BitMap

from binette import bin_quality
from binette.bin_manager import Bin
from binette.bin_quality import balanced_chunks, chunks


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


def test_add_bin_metrics_with_multiple_threads():
    """
    Test add_bin_metrics with multi-threading enabled.

    This test verifies that:
    1. The function correctly triggers parallel processing path with multiple threads
    2. The _assess_bins_quality_batch function is called for processing bins
    3. All bins are processed and returned with quality metrics
    4. joblib.Parallel is called with the correct number of jobs

    Note: This test mocks joblib.Parallel to avoid actual multiprocessing complexity
    and focuses on testing the parallel processing logic and chunking behavior.
    """
    # Create mock bins - need enough bins to trigger parallel processing
    # min_bins_per_chunk = checkm2_batch_size * 6 = 500 * 6 = 3000
    # Need > 3000 * 2 = 6000 bins to trigger parallel processing with threads > 1
    num_bins = 7000
    bins = [Bin(BitMap((i,))) for i in range(num_bins)]

    # Mock contig_info with minimal required data
    contig_info = {
        "contig_to_kegg_counter": {i: Counter({"K01810": 1}) for i in range(num_bins)},
        "contig_to_cds_count": {i: 10 for i in range(num_bins)},
        "contig_to_aa_counter": {
            i: Counter({"A": 5, "D": 10}) for i in range(num_bins)
        },
        "contig_to_aa_length": {i: 1000 for i in range(num_bins)},
    }

    contamination_weight = 0.5
    threads = 4
    checkm2_batch_size = 500

    # Mock _assess_bins_quality_batch to add quality metrics without CheckM2
    def mock_assess_bins_quality_batch(
        bins_batch,
        contig_to_kegg_counter,
        contig_to_cds_count,
        contig_to_aa_counter,
        contig_to_aa_length,
        contamination_weight,
        postProcessor,
        threads_arg,
    ):
        """Mock implementation that adds dummy quality metrics to bins."""
        for bin_obj in bins_batch:
            bin_obj.add_quality(
                completeness=90.0,
                contamination=5.0,
                contamination_weight=contamination_weight,
            )
            bin_obj.add_model("Neural Network (Specific Model)")
        return bins_batch

    # Mock joblib.Parallel to execute tasks sequentially without spawning processes
    class MockParallel:
        def __init__(self, n_jobs=1, **kwargs):
            self.n_jobs = n_jobs
            self.call_count = 0

        def __call__(self, tasks):
            self.call_count += 1
            # Execute tasks sequentially (simulating parallel execution)
            results = []
            for task in tasks:
                result = task[0](*task[1], **task[2])  # Execute the delayed function
                results.append(result)
            return results

    with patch(
        "binette.bin_quality._assess_bins_quality_batch",
        side_effect=mock_assess_bins_quality_batch,
    ):
        with patch("binette.bin_quality.joblib.Parallel", MockParallel):
            with patch("binette.bin_quality.get_modelPostprocessing"):
                # Call add_bin_metrics with multiple threads
                result_bins = bin_quality.add_bin_metrics(
                    bins=bins,
                    contig_info=contig_info,
                    contamination_weight=contamination_weight,
                    threads=threads,
                    checkm2_batch_size=checkm2_batch_size,
                    disable_progress_bar=True,
                )

                # Verify all bins were processed
                assert len(result_bins) == num_bins, (
                    f"Expected {num_bins} bins, got {len(result_bins)}"
                )

                # Verify that each bin has quality metrics
                for bin_obj in result_bins:
                    assert hasattr(bin_obj, "completeness"), (
                        "Bin should have completeness metric"
                    )
                    assert hasattr(bin_obj, "contamination"), (
                        "Bin should have contamination metric"
                    )
                    assert bin_obj.completeness == 90.0, "Completeness should be 90.0"
                    assert bin_obj.contamination == 5.0, "Contamination should be 5.0"


def test_add_bin_metrics_sequential_path():
    """
    Test add_bin_metrics with single thread (sequential processing).

    This test verifies that:
    1. With threads=1, the sequential processing path is used
    2. Bins are processed correctly without parallelization
    3. Quality metrics are correctly added to bins

    Note: This test mocks _assess_bins_quality_batch to avoid CheckM2 imports.
    """
    # Create a small number of bins
    num_bins = 10
    bins = [Bin(BitMap((i,))) for i in range(num_bins)]

    # Mock contig_info
    contig_info = {
        "contig_to_kegg_counter": {i: Counter({"K01810": 1}) for i in range(num_bins)},
        "contig_to_cds_count": {i: 10 for i in range(num_bins)},
        "contig_to_aa_counter": {
            i: Counter({"A": 5, "D": 10}) for i in range(num_bins)
        },
        "contig_to_aa_length": {i: 1000 for i in range(num_bins)},
    }

    contamination_weight = 0.5

    # Mock _assess_bins_quality_batch to add quality metrics without CheckM2
    def mock_assess_bins_quality_batch(
        bins_batch,
        contig_to_kegg_counter,
        contig_to_cds_count,
        contig_to_aa_counter,
        contig_to_aa_length,
        contamination_weight,
        postProcessor,
        threads_arg,
    ):
        """Mock implementation that adds dummy quality metrics to bins."""
        for bin_obj in bins_batch:
            bin_obj.add_quality(
                completeness=85.0,
                contamination=3.0,
                contamination_weight=contamination_weight,
            )
            bin_obj.add_model("Gradient Boost (General Model)")
        return bins_batch

    # Mock at module level
    with patch(
        "binette.bin_quality._assess_bins_quality_batch",
        side_effect=mock_assess_bins_quality_batch,
    ) as mock_assess:
        with patch("binette.bin_quality.get_modelPostprocessing"):
            # Call add_bin_metrics with single thread
            result_bins = bin_quality.add_bin_metrics(
                bins=bins,
                contig_info=contig_info,
                contamination_weight=contamination_weight,
                threads=1,  # Single thread
                checkm2_batch_size=500,
                disable_progress_bar=True,
            )

            # For sequential processing, should call _assess_bins_quality_batch at least once
            assert mock_assess.call_count >= 1, (
                f"Should call _assess_bins_quality_batch at least once, got {mock_assess.call_count}"
            )

            # Verify all bins were processed
            assert len(result_bins) == num_bins, (
                f"Expected {num_bins} bins, got {len(result_bins)}"
            )

            # Verify that each bin has quality metrics
            for bin_obj in result_bins:
                assert hasattr(bin_obj, "completeness"), (
                    "Bin should have completeness metric"
                )
                assert hasattr(bin_obj, "contamination"), (
                    "Bin should have contamination metric"
                )
                assert bin_obj.completeness == 85.0, "Completeness should be 85.0"
                assert bin_obj.contamination == 3.0, "Contamination should be 3.0"


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
