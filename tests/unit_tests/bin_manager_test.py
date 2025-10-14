"""
Unit tests for binette.

"""

import logging
from pathlib import Path

import networkx as nx
import pytest
from pyroaring import BitMap

from binette import bin_manager, bin_quality


def test_get_all_possible_combinations():
    input_list = ["2", "3", "4"]
    expected_list = [("2", "3"), ("2", "4"), ("3", "4"), ("2", "3", "4")]

    assert list(bin_manager.get_all_possible_combinations(input_list)) == expected_list


def test_bin_eq_true():
    bin1 = bin_manager.Bin(contigs=BitMap({1, 2}), origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs=BitMap({1, 2}), origin="test2", name="binA")

    assert bin1 == bin2


def test_bin_eq_false():
    bin1 = bin_manager.Bin(contigs=BitMap({1, 2}), origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs=BitMap({1, 2, 3}), origin="test2", name="binA")

    assert bin1 != bin2


def test_in_for_bin_list():
    bin1 = bin_manager.Bin(contigs=BitMap({1, 2, 3}), origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs=BitMap({1, 2, 33}), origin="test2", name="binA")
    bin3 = bin_manager.Bin(contigs=BitMap({4, 5}), origin="test2", name="binA")

    bins = [bin1, bin2]

    assert bin1 in bins
    assert bin2 in bins
    assert bin3 not in bins


def test_add_length_positive_integer():
    bin_obj = bin_manager.Bin(contigs=BitMap({1, 2, 3}), origin="test1", name="bin1")
    length = 100
    bin_obj.add_length(length)
    assert bin_obj.length == length


def test_add_length_negative_integer():
    bin_obj = bin_manager.Bin(contigs=BitMap({1, 2, 3}), origin="test1", name="bin1")
    with pytest.raises(ValueError):
        length = -50
        bin_obj.add_length(length)


def test_add_n50_positive_integer():
    bin_obj = bin_manager.Bin(contigs=BitMap({1, 2, 3}), origin="test1", name="bin1")
    n50 = 100
    bin_obj.add_N50(n50)
    assert bin_obj.N50 == n50


def test_add_n50_negative_integer():
    bin_obj = bin_manager.Bin(contigs=BitMap({1, 2, 3}), origin="test1", name="bin1")
    with pytest.raises(ValueError):
        n50 = -50
        bin_obj.add_N50(n50)


def test_add_quality():
    completeness = 10
    contamination = 6
    contamination_weight = 2

    bin_obj = bin_manager.Bin(contigs=BitMap({1, 2, 3}), origin="test1", name="bin1")

    bin_obj.add_quality(completeness, contamination, contamination_weight)

    assert bin_obj.completeness == completeness
    assert bin_obj.contamination == contamination

    assert bin_obj.score == completeness - contamination * contamination_weight


def test_add_model():
    bin_obj = bin_manager.Bin(contigs=BitMap({1, 2, 3}), origin="test1", name="bin1")

    bin_obj.add_model("Neural Network (Specific Model)")

    assert bin_obj.checkm2_model == "Neural Network (Specific Model)"


def test_add_model_error():
    bin_obj = bin_manager.Bin(contigs=BitMap({1, 2, 3}), origin="test1", name="bin1")

    with pytest.raises(ValueError):
        bin_obj.add_model("Not a valid model name")


def test_is_high_quality():
    completeness = 90
    contamination = 1
    contamination_weight = 2

    bin_obj = bin_manager.Bin(contigs=BitMap({1, 2, 3}), origin="test1", name="bin1")

    bin_obj.add_quality(completeness, contamination, contamination_weight)

    assert bin_obj.is_high_quality(min_completeness=80, max_contamination=5) is True


def test_is_high_quality_no_quality():
    bin_obj = bin_manager.Bin(contigs=BitMap({1, 2, 3}), origin="test1", name="bin1")

    with pytest.raises(ValueError):
        bin_obj.is_high_quality(min_completeness=80, max_contamination=5)


def test_multiple_bins_intersection():
    bin1 = bin_manager.Bin(contigs=BitMap({1, 2, 3, 987}), origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs=BitMap({1, 2, 33}), origin="test2", name="binA")
    bin3 = bin_manager.Bin(contigs=BitMap({1, 2, 33}), origin="test2", name="binA")

    bin_intersection = bin1.contig_intersection(bin2, bin3)

    assert bin_intersection == BitMap({1, 2})


def test_bin_overlap_true():
    bin1 = bin_manager.Bin(contigs=BitMap({1, 2, 3}), origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs=BitMap({1, 2, 33}), origin="test2", name="binA")

    assert bin1.overlaps_with(bin2)

    assert bin2.overlaps_with(bin1)


def test_bin_overlap_false():
    bin1 = bin_manager.Bin(contigs=BitMap({13, 21, 37}), origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs=BitMap({1, 2, 33}), origin="test2", name="binA")

    assert not bin1.overlaps_with(bin2)


def test_bin_union():
    bin1 = bin_manager.Bin(contigs=BitMap({13, 21}), origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs=BitMap({1, 2, 33}), origin="test2", name="binA")

    expected_union_bin_composition = BitMap({13, 21, 1, 2, 33})

    union_bin = bin1.contig_union(bin2)

    assert union_bin == expected_union_bin_composition


def test_bin_union2():
    # Create some example bins
    bin1 = bin_manager.Bin(BitMap({1, 2}), "origin1", "bin1")
    bin2 = bin_manager.Bin(BitMap({2, 3}), "origin2", "bin2")
    bin3 = bin_manager.Bin(BitMap({4, 5}), "origin3", "bin3")

    # Perform union operation
    union_bin = bin1.contig_union(bin2, bin3)

    # Check the result
    expected_contigs = BitMap({1, 2, 3, 4, 5})

    assert union_bin == expected_contigs


def test_no_bitmap():
    with pytest.raises(TypeError):
        bin_manager.Bin({1, 2}, "origin1", "bin1")


def test_bin_difference():
    bin1 = bin_manager.Bin(contigs=BitMap({5, 1, 6, 7, 8}), origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs=BitMap({3, 6, 7}), origin="test2", name="bin2")
    bin3 = bin_manager.Bin(contigs=BitMap({1, 2, 3, 4, 6}), origin="test2", name="bin3")

    diff_bin1_2_3 = BitMap({6})
    diff_bin1_2 = BitMap({6, 7})

    assert bin1.contig_intersection(bin2, bin3) == diff_bin1_2_3
    assert bin1.contig_intersection(bin2) == diff_bin1_2


def test_bin_intersection():
    bin1 = bin_manager.Bin(contigs=BitMap({5, 1, 6, 7, 8}), origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs=BitMap({3, 6, 7}), origin="test2", name="bin2")
    bin3 = bin_manager.Bin(contigs=BitMap({1, 2, 3, 4, 7}), origin="test2", name="bin3")

    inter_bin123 = BitMap({7})
    inter_bin1_2 = BitMap({7, 6})

    assert bin1.contig_intersection(bin2, bin3) == inter_bin123
    assert bin1.contig_intersection(bin2) == inter_bin1_2


def test_select_best_bins_simple():
    b1 = bin_manager.Bin(contigs=BitMap({1, 2}), origin="", name="")
    b2 = bin_manager.Bin(contigs=BitMap({2}), origin="", name="")
    b3 = bin_manager.Bin(contigs=BitMap({3}), origin="", name="")

    b1.score = 90
    b2.score = 80
    b3.score = 70

    b1.completeness = 96
    b2.completeness = 96
    b3.completeness = 96

    b1.contamination = 0.5
    b2.contamination = 5
    b3.contamination = 0.5

    b1.N50 = 100
    b2.N50 = 100
    b3.N50 = 100

    assert bin_manager.select_best_bins(
        {1: b1, 2: b2, 3: b3}, min_completeness=40, max_contamination=10
    ) == [b1, b3]


def test_select_best_bins_with_same_score():
    b1 = bin_manager.Bin(contigs=BitMap({1, 2}), origin="", name="")
    b2 = bin_manager.Bin(contigs=BitMap({2}), origin="", name="")
    b3 = bin_manager.Bin(contigs=BitMap({3}), origin="", name="")

    b1.score = 90
    b2.score = 90
    b3.score = 70

    b1.N50 = 100
    b2.N50 = 101  # selection is then based on N50
    b3.N50 = 100

    b1.completeness = 90
    b2.completeness = 90
    b3.completeness = 70

    b1.contamination = 0
    b2.contamination = 0
    b3.contamination = 0

    assert bin_manager.select_best_bins(
        {1: b1, 2: b2, 3: b3}, max_contamination=0, min_completeness=40
    ) == [b2, b3]


def test_select_best_bins_with_equality_based_on_original():
    b1 = bin_manager.Bin(contigs=BitMap({1, 2}), origin="", name="", is_original=False)
    b2 = bin_manager.Bin(contigs=BitMap({2}), origin="", name="", is_original=True)
    b3 = bin_manager.Bin(contigs=BitMap({3}), origin="", name="", is_original=False)

    b1.score = 90
    b2.score = 90
    b3.score = 70

    b1.N50 = 100
    b2.N50 = 100
    b3.N50 = 100

    b1.completeness = 90
    b2.completeness = 90
    b3.completeness = 70

    b1.contamination = 0
    b2.contamination = 0
    b3.contamination = 0
    # when score and n50 is the same, selection is made on the smallest key.
    # bin created first have a smaller id. so b1 should selected
    assert bin_manager.select_best_bins(
        {1: b1, 2: b2, 3: b3}, max_contamination=0, min_completeness=40
    ) == [b2, b3]


def test_select_best_bins_with_equality():
    b1 = bin_manager.Bin(contigs=BitMap({1, 2}), origin="", name="", is_original=False)
    b2 = bin_manager.Bin(contigs=BitMap({2}), origin="", name="", is_original=False)
    b3 = bin_manager.Bin(contigs=BitMap({3}), origin="", name="", is_original=False)

    b1.score = 90
    b2.score = 90
    b3.score = 70

    b1.N50 = 100
    b2.N50 = 100
    b3.N50 = 100
    b1.completeness = 90
    b2.completeness = 90
    b3.completeness = 70

    b1.contamination = 0
    b2.contamination = 0
    b3.contamination = 0
    # when score, n50 and is_original is the same, selection is made on the smallest key.
    # bin created first have a smaller id. so b1 should selected
    assert bin_manager.select_best_bins(
        {1: b1, 2: b2, 3: b3}, max_contamination=0, min_completeness=40
    ) == [b1, b3]


# The function should create intersection bins when there are overlapping contigs between bins.
def test_intersection_bins_created():
    set1 = [
        bin_manager.Bin(contigs=BitMap({1, 2}), origin="A", name="bin1"),
        bin_manager.Bin(contigs=BitMap({3, 4}), origin="A", name="bin2"),
        bin_manager.Bin(contigs=BitMap({5}), origin="A", name="bin2"),
    ]
    # need to defined completeness and conta
    # because when too low the bin is not used in all operation
    for b in set1:
        b.completeness = 100
        b.contamination = 0
        b.length = 7000

    binA = bin_manager.Bin(contigs=BitMap({1, 3}), origin="B", name="binA")
    binA.contamination = 0
    binA.completeness = 100
    binA.length = 7000
    set2 = [
        binA,
    ]
    input_bins = set1 + set2

    key_to_bins = {b.contigs_key: b for b in input_bins}

    contig_lengths = bin_quality.prepare_contig_sizes(
        {1: 5000, 2: 3000, 3: 4000, 4: 2000, 5: 1000}
    )

    intermediate_bins_result = bin_manager.create_intermediate_bins(
        key_to_bins,
        contig_lengths=contig_lengths,
        min_len=0,
        max_len=10_000_000,
        min_comp=0,
        max_conta=100,
    )

    expected_bin_compositions = [
        BitMap({1, 2, 3}),
        BitMap({2}),
        BitMap({1}),
        BitMap({1, 4, 3}),
        BitMap({4}),
        BitMap({3}),
    ]

    for b in intermediate_bins_result.values():
        assert b.contigs in expected_bin_compositions
    assert len(intermediate_bins_result) == len(expected_bin_compositions)


def test_from_bins_to_bin_graph():
    bin1 = bin_manager.Bin(contigs=BitMap({1, 2}), origin="A", name="bin1")
    bin2 = bin_manager.Bin(contigs=BitMap({3, 4}), origin="A", name="bin2")
    bin3 = bin_manager.Bin(contigs=BitMap({5}), origin="A", name="bin3")

    set1 = [bin1, bin2, bin3]

    binA = bin_manager.Bin(contigs=BitMap({1, 3}), origin="B", name="binA")

    set2 = [binA]

    result_graph = bin_manager.from_bins_to_bin_graph(set1 + set2)

    assert result_graph.number_of_edges() == 2
    # bin3 is not connected to any bin so it is not in the graph
    assert result_graph.number_of_nodes() == 3

    assert set(result_graph.nodes) == {b.contigs_key for b in [binA, bin1, bin2]}


@pytest.fixture
def simple_bin_graph():
    bin1 = bin_manager.Bin(contigs=BitMap({1, 2, 3}), origin="A", name="bin1")
    bin2 = bin_manager.Bin(contigs=BitMap({1, 2, 4}), origin="B", name="bin2")

    for b in [bin1, bin2]:
        b.completeness = 100
        b.contamination = 0

    G = nx.Graph()
    G.add_edge(bin1, bin2)

    return G


def test_get_bins_from_contig2bin_table(tmp_path):
    # Create a temporary file (contig-to-bin table) for testing
    test_table_content = [
        "# Sample contig-to-bin table",
        "contig1\tbin1",
        "contig2\tbin1",
        "contig3\tbin2",
    ]
    test_table_path = tmp_path / "test_contig2bin_table.txt"
    test_table_path.write_text("\n".join(test_table_content))

    # Define set name for the bins
    set_name = "TestSet"

    # Call the function to generate Bin objects
    result_bins = bin_manager.get_bins_from_contig2bin_table(
        str(test_table_path), set_name
    )

    # Validate the result
    assert len(result_bins) == 2  # Check if the correct number of bins are created

    # Define expected bins based on the test table content
    expected_bins = [
        {"contigs": {"contig1", "contig2"}, "set_name": set_name, "bin_name": "bin1"},
        {"contigs": {"contig3"}, "set_name": set_name, "bin_name": "bin2"},
    ]

    # Compare expected bins with the result
    assert all(expected_bin in result_bins for expected_bin in expected_bins)
    assert all(result_bin in expected_bins for result_bin in result_bins)


def test_parse_contig2bin_tables(tmp_path):
    # Create temporary contig-to-bin tables for testing
    test_tables = {
        "set1": [
            "# Sample contig-to-bin table for bin1",
            "contig1\tbin1",
            "contig2\tbin1",
            "contig3\tbin2",
        ],
        "set2": [
            "# Sample contig-to-bin table for bin2",
            "contig3\tbinA",
            "contig4\tbinA",
        ],
    }

    # Create temporary files for contig-to-bin tables
    for name, content in test_tables.items():
        table_path = tmp_path / f"test_{name}_contig2bin_table.txt"
        table_path.write_text("\n".join(content))

    # Call the function to parse contig-to-bin tables
    result_bin_dict = bin_manager.parse_contig2bin_tables(
        {
            name: str(tmp_path / f"test_{name}_contig2bin_table.txt")
            for name in test_tables
        }
    )

    # Validate the result
    assert len(result_bin_dict) == len(
        test_tables
    )  # Check if the number of bins matches the number of tables

    # Define expected Bin objects based on the test tables
    expected_bins = {
        "set1": [
            {"contigs": {"contig1", "contig2"}, "set_name": "set1", "bin_name": "bin1"},
            {"contigs": {"contig3"}, "set_name": "set1", "bin_name": "bin2"},
        ],
        "set2": [
            {"contigs": {"contig3", "contig4"}, "set_name": "set2", "bin_name": "binA"},
        ],
    }

    # Compare expected bins with the result
    for name, expected in expected_bins.items():
        assert name in result_bin_dict
        assert len(result_bin_dict[name]) == len(expected)
        for result_bin in result_bin_dict[name]:
            assert result_bin in expected


@pytest.fixture
def create_temp_bin_files(tmpdir):
    # Create temporary bin files
    bin_dir = tmpdir.mkdir("bins")
    bin1 = bin_dir.join("bin1.fasta")
    bin1.write(">contig1\nATGC\n>contig2\nGCTA")

    bin2 = bin_dir.join("bin2.fasta")
    bin2.write(">contig3\nTTAG\n>contig4\nCGAT")

    return bin_dir


@pytest.fixture
def create_temp_bin_directories(tmpdir, create_temp_bin_files):
    # Create temporary bin directories
    bin_dir1 = tmpdir.mkdir("set1")
    bin1 = bin_dir1.join("bin1.fasta")
    bin1.write(">contig1\nATGC\n>contig2\nGCTA")

    bin2 = bin_dir1.join("bin2.fasta")
    bin2.write(">contig3\nTTAG\n>contig4\nCGAT")

    bin_dir2 = tmpdir.mkdir("set2")
    bin2 = bin_dir2.join("binA.fasta")
    bin2.write(">contig3\nTTAG\n>contig4\nCGAT\n>contig5\nCGGC")

    return {"set1": Path(bin_dir1), "set2": Path(bin_dir2)}


def test_get_bins_from_directory(create_temp_bin_files):
    bin_dir = create_temp_bin_files
    set_name = "TestSet"

    bins = bin_manager.get_bins_from_directory(
        Path(bin_dir), set_name, fasta_extensions={".fasta"}
    )

    assert len(bins) == 2  # Ensure that the correct number of Bin objects is returned

    # Check if the Bin objects are created with the correct contigs, set name, and bin names
    assert isinstance(bins[0], dict)
    assert isinstance(bins[1], dict)
    assert bins[1]["contigs"] in [{"contig1", "contig2"}, {"contig3", "contig4"}]
    assert bins[0]["contigs"] in [{"contig1", "contig2"}, {"contig3", "contig4"}]
    assert bins[0]["set_name"] == set_name
    assert bins[1]["set_name"] == set_name
    assert bins[1]["bin_name"] in ["bin2", "bin1"]
    assert bins[0]["bin_name"] in ["bin2", "bin1"]


def test_get_bins_from_directory_no_files(tmpdir):
    bin_dir = Path(tmpdir.mkdir("empty_bins"))
    set_name = "EmptySet"

    bins = bin_manager.get_bins_from_directory(
        bin_dir, set_name, fasta_extensions={".fasta"}
    )

    assert (
        len(bins) == 0
    )  # Ensure that no Bin objects are returned for an empty directory


def test_get_bins_from_directory_no_wrong_extensions(create_temp_bin_files):
    bin_dir = Path(create_temp_bin_files)
    set_name = "TestSet"

    bins = bin_manager.get_bins_from_directory(
        bin_dir, set_name, fasta_extensions={".fna"}
    )

    assert (
        len(bins) == 0
    )  # Ensure that no Bin objects are returned for an empty directory


def test_parse_bin_directories(create_temp_bin_directories):
    set_name_to_bin_dir = create_temp_bin_directories

    bins = bin_manager.parse_bin_directories(
        set_name_to_bin_dir, fasta_extensions={".fasta"}
    )

    assert len(bins) == 2  # Ensure that the correct number of bin directories is parsed

    # Check if the Bin objects are created with the correct contigs, set name, and bin names
    assert isinstance(list(bins["set1"])[0], dict)
    assert isinstance(list(bins["set2"])[0], dict)

    assert len(bins["set2"]) == 1
    assert len(bins["set1"]) == 2


def test_get_contigs_in_bin_sets(caplog):
    """
    Test the get_contigs_in_bin_sets function for correct behavior.

    :param mock_bins: The mock_bins fixture providing test bin data.
    :param caplog: The pytest caplog fixture to capture logging output.
    """

    bin_set_name_to_bins = {
        "set1": [
            {"contigs": {"1", "4"}, "set_name": "set1", "bin_name": "binA"},
            {
                "contigs": {"1"},
                "set_name": "set1",
                "bin_name": "bin1",
            },
        ],
        "set2": [{"contigs": {"3", "4", "5"}, "set_name": "set2", "bin_name": "binB"}],
    }

    # Test the function with valid data
    # warning because set1 has duplicated contig "1"
    with caplog.at_level(logging.WARNING):
        result = bin_manager.get_contigs_in_bin_sets(bin_set_name_to_bins)

    # Expected unique contigs
    expected_contigs = {"1", "3", "4", "5"}

    # Check if the result matches expected contigs
    assert set(result) == expected_contigs, "The returned set of contigs is incorrect."
