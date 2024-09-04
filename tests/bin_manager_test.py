"""
Unit tests for binette.

"""

import pytest

from binette import bin_manager
import networkx as nx

import logging
from pathlib import Path

def test_get_all_possible_combinations():
    input_list = ["2", "3", "4"]
    expected_list = [("2", "3"), ("2", "4"), ("3", "4"), ("2", "3", "4")]

    assert list(bin_manager.get_all_possible_combinations(input_list)) == expected_list


@pytest.fixture
def example_bin_set1():
    bin1 = bin_manager.Bin(contigs={"1", "2"}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={"3", "4"}, origin="test1", name="bin2")
    bin3 = bin_manager.Bin(contigs={"5"}, origin="test1", name="bin2")
    return {bin1, bin2, bin3}

@pytest.fixture
def example_bin_set2():
    bin1 = bin_manager.Bin(contigs={"1", "2", "3"}, origin="test2", name="binA")
    bin2 = bin_manager.Bin(contigs={"4", "5"}, origin="test2", name="binB")
    return {bin1, bin2}

def test_bin_eq_true():
    bin1 = bin_manager.Bin(contigs={"1", "2", "e"}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2"}, origin="test2", name="binA")

    assert bin1 == bin2


def test_bin_eq_false():
    bin1 = bin_manager.Bin(contigs={"1", "2", "e"}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin="test2", name="binA")

    assert bin1 != bin2


def test_in_for_bin_list():
    bin1 = bin_manager.Bin(contigs={"1", "2", "e"}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin="test2", name="binA")
    bin3 = bin_manager.Bin(contigs={"4", "R"}, origin="test2", name="binA")

    bins = [bin1, bin2]

    assert bin1 in bins
    assert bin2 in bins
    assert bin3 not in bins

def test_add_length_positive_integer():
    bin_obj = bin_manager.Bin(contigs={"1", "2", "e"}, origin="test1", name="bin1")
    length = 100
    bin_obj.add_length(length)
    assert bin_obj.length == length

def test_add_length_negative_integer():
    bin_obj = bin_manager.Bin(contigs={"1", "2", "e"}, origin="test1", name="bin1") 
    with pytest.raises(ValueError):
        length = -50
        bin_obj.add_length(length)

def test_add_n50_positive_integer():
    bin_obj = bin_manager.Bin(contigs={"1", "2", "e"}, origin="test1", name="bin1")
    n50 = 100
    bin_obj.add_N50(n50)
    assert bin_obj.N50 == n50

def test_add_n50_negative_integer():
    bin_obj = bin_manager.Bin(contigs={"1", "2", "e"}, origin="test1", name="bin1") 
    with pytest.raises(ValueError):
        n50 = -50
        bin_obj.add_N50(n50)

def test_add_quality():

    completeness = 10
    contamination = 6
    contamination_weight = 2

    bin_obj = bin_manager.Bin(contigs={"1", "2", "e"}, origin="test1", name="bin1") 

    bin_obj.add_quality(completeness, contamination, contamination_weight)
    
    assert bin_obj.completeness == completeness 
    assert bin_obj.contamination == contamination 

    assert bin_obj.score == completeness - contamination * contamination_weight



    


# def test_two_bin_intersection():
#     bin1 = bin_manager.Bin(contigs={"1", "2", "e", "987"}, origin="test1", name="bin1")
#     bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin="test2", name="binA")

#     bin_intersection = bin1 & bin2

#     assert bin_intersection == bin_manager.Bin({"1", "2", "e"}, "", "")


def test_multiple_bins_intersection():
    bin1 = bin_manager.Bin(contigs={"1", "2", "e", "987"}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin="test2", name="binA")
    bin3 = bin_manager.Bin(contigs={"1", "123", "2", "33"}, origin="test2", name="binA")

    bin_intersection = bin1.intersection(bin2, bin3)

    assert bin_intersection == bin_manager.Bin({"1", "2"}, "", "")


def test_bin_set_has_no_duplicate(example_bin_set1):
    set1_twice = list(example_bin_set1) + list(example_bin_set1)

    assert example_bin_set1 == set(set1_twice)


def test_bin_overlap_true():
    bin1 = bin_manager.Bin(contigs={"1", "2", "e"}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin="test2", name="binA")

    assert bin1.overlaps_with(bin2)

    assert bin2.overlaps_with(bin1)


def test_bin_overlap_false():
    bin1 = bin_manager.Bin(contigs={"13", "21", "ef"}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin="test2", name="binA")

    assert not bin1.overlaps_with(bin2)


def test_bin_union():
    bin1 = bin_manager.Bin(contigs={"13", "21"}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin="test2", name="binA")

    expected_union_bin = bin_manager.Bin(contigs={"13", "21", "1", "e", "2", "33"}, origin="", name="")
    union_bin = bin1.union(bin2)

    assert union_bin == expected_union_bin
    assert union_bin.name == f"{bin1.id} | {bin2.id}"


def test_bin_union2():
    # Create some example bins
    bin1 = bin_manager.Bin({'contig1', 'contig2'}, 'origin1', 'bin1')
    bin2 = bin_manager.Bin({'contig2', 'contig3'}, 'origin2', 'bin2')
    bin3 = bin_manager.Bin({'contig4', 'contig5'}, 'origin3', 'bin3')

    # Perform union operation
    union_bin = bin1.union(bin2, bin3)

    # Check the result
    expected_contigs = {'contig1', 'contig2', 'contig3', 'contig4', 'contig5'}
    expected_origin = {'union'}

    assert union_bin.contigs == expected_contigs
    assert union_bin.origin == expected_origin


def test_bin_difference():
    bin1 = bin_manager.Bin(contigs={5, 1, 6, 7, 8}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={3, 6, 7}, origin="test2", name="bin2")
    bin3 = bin_manager.Bin(contigs={1, 2, 3, 4}, origin="test2", name="bin3")

    diff_bin1_23 = bin_manager.Bin(contigs={5, 8}, origin="", name="")
    diff_bin1_2 = bin_manager.Bin(contigs={5, 1, 8}, origin="", name="")

    assert bin1.difference(bin2, bin3) == diff_bin1_23
    assert bin1.difference(bin2) == diff_bin1_2
    assert bin1.difference(bin2, bin3).name == f"{bin1.id} - {bin2.id} - {bin3.id}"


def test_bin_intersection():
    bin1 = bin_manager.Bin(contigs={5, 1, 6, 7, 8}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={3, 6, 7}, origin="test2", name="bin2")
    bin3 = bin_manager.Bin(contigs={1, 2, 3, 4, 7}, origin="test2", name="bin3")

    inter_bin123 = bin_manager.Bin(contigs={7}, origin="", name="")
    iner_bin1_2 = bin_manager.Bin(contigs={7, 6}, origin="", name="")

    assert bin1.intersection(bin2, bin3) == inter_bin123
    assert bin1.intersection(bin2) == iner_bin1_2
    assert bin1.intersection(bin2, bin3).name == f"{bin1.id} & {bin2.id} & {bin3.id}"


def test_select_best_bins_simple():

    b1 = bin_manager.Bin(contigs={1, 2}, origin="", name="")
    b2 = bin_manager.Bin(contigs={2}, origin="", name="")
    b3 = bin_manager.Bin(contigs={3}, origin="", name="")

    b1.score = 90
    b2.score = 80
    b3.score = 70

    b1.N50 = 100
    b2.N50 = 100
    b3.N50 = 100

    assert bin_manager.select_best_bins({b1, b2, b3}) == [b1, b3]


def test_select_best_bins_with_same_score():

    b1 = bin_manager.Bin(contigs={1, 2}, origin="", name="")
    b2 = bin_manager.Bin(contigs={2}, origin="", name="")
    b3 = bin_manager.Bin(contigs={3}, origin="", name="")

    b1.score = 90
    b2.score = 90
    b3.score = 70

    b1.N50 = 100
    b2.N50 = 101
    b3.N50 = 100

    assert bin_manager.select_best_bins({b1, b2, b3}) == [b2, b3]


def test_select_best_bins_with_equality():

    b1 = bin_manager.Bin(contigs={1, 2}, origin="", name="")
    b2 = bin_manager.Bin(contigs={2}, origin="", name="")
    b3 = bin_manager.Bin(contigs={3}, origin="", name="")

    b1.score = 90
    b2.score = 90
    b3.score = 70

    b1.N50 = 100
    b2.N50 = 100
    b3.N50 = 100

    # when score and n50 is the same, selection is made on the smallest id.
    # bin created first have a smaller id. so b1 should selected
    assert bin_manager.select_best_bins({b1, b2, b3}) == [b1, b3]


# The function should create intersection bins when there are overlapping contigs between bins.
def test_intersection_bins_created():
    set1 = [
        bin_manager.Bin(contigs={"1", "2"}, origin="A", name="bin1"),
        bin_manager.Bin(contigs={"3", "4"}, origin="A", name="bin2"),
        bin_manager.Bin(contigs={"5"}, origin="A", name="bin2"),
    ]
    # need to defined completeness and conta
    # because when too low the bin is not used in all operation
    for b in set1:
        b.completeness = 100
        b.contamination = 0

    binA = bin_manager.Bin(contigs={"1", "3"}, origin="B", name="binA")
    binA.contamination = 0
    binA.completeness = 100
    set2 = [
        binA,
    ]
    bin_set_name_to_bins = {'set1': set1,
                            "set2":set2}

    intermediate_bins_result = bin_manager.create_intermediate_bins(bin_set_name_to_bins)

    expected_intermediate_bins = {bin_manager.Bin(contigs={"1", "2", "3"}, origin="bin1 | binA ", name="NA"),
                                  bin_manager.Bin(contigs={"2"}, origin="bin1 - binA ", name="NA"),
                                  bin_manager.Bin(contigs={"1"}, origin="bin1 & binA ", name="NA"),
                                  bin_manager.Bin(contigs={"1", "4", "3"}, origin="bin2 | binA ", name="NA"),
                                  bin_manager.Bin(contigs={"4"}, origin="bin2 - binA ", name="NA"), # binA -bin1 is equal to bin1 & binA
                                  bin_manager.Bin(contigs={"3"}, origin="bin1 & binA ", name="NA"),
                                  }

    assert intermediate_bins_result == expected_intermediate_bins


# Renames contigs in bins based on provided mapping.
def test_renames_contigs(example_bin_set1):
    
    bin_set = [
        bin_manager.Bin(contigs={"c1", "c2"}, origin="A", name="bin1"),
        bin_manager.Bin(contigs={"c3", "c4"}, origin="A", name="bin2")
    ]
        
    contig_to_index = {'c1': 1, 'c2': 2, 'c3': 3, 'c4': 4, "c5":5}

    # Act
    bin_manager.rename_bin_contigs(bin_set, contig_to_index)

    # Assert
    assert bin_set[0].contigs == {1, 2}
    assert bin_set[0].hash == hash(str(sorted({1, 2})))
    assert bin_set[1].contigs == {3, 4}
    assert bin_set[1].hash == hash(str(sorted({3, 4})))

def test_get_contigs_in_bins():
    bin_set = [
        bin_manager.Bin(contigs={"c1", "c2"}, origin="A", name="bin1"),
        bin_manager.Bin(contigs={"c3", "c4"}, origin="A", name="bin2"),
        bin_manager.Bin(contigs={"c3", "c18"}, origin="A", name="bin2")
    ]

    contigs = bin_manager.get_contigs_in_bins(bin_set)

    assert set(contigs) == {"c1", "c2", "c3", "c4", "c18"}


def test_dereplicate_bin_sets():
    b1 = bin_manager.Bin(contigs={"c1", "c2"}, origin="A", name="bin1")
    b2 = bin_manager.Bin(contigs={"c3", "c4"}, origin="A", name="bin2")
    b3 = bin_manager.Bin(contigs={"c3", "c18"}, origin="A", name="bin2")

    bdup  = bin_manager.Bin(contigs={"c3", "c18"}, origin="D", name="C")

    derep_bins_result =  bin_manager.dereplicate_bin_sets([[b2,b3 ],[b1, bdup]])


    assert derep_bins_result == {b1, b2, b3}



def test_from_bin_sets_to_bin_graph():


    bin1 = bin_manager.Bin(contigs={"1", "2"}, origin="A", name="bin1")
    bin2 = bin_manager.Bin(contigs={"3", "4"}, origin="A", name="bin2")
    bin3 = bin_manager.Bin(contigs={"5"}, origin="A", name="bin3")

    set1 = [bin1,bin2,bin3]

    binA = bin_manager.Bin(contigs={"1", "3"}, origin="B", name="binA")
            
    set2 = [binA]

    result_graph = bin_manager.from_bin_sets_to_bin_graph({"B":set2, "A":set1}) 

    assert result_graph.number_of_edges() == 2
    # bin3 is not connected to any bin so it is not in the graph
    assert result_graph.number_of_nodes() == 3

    assert set(result_graph.nodes) == {binA, bin1, bin2}

@pytest.fixture
def simple_bin_graph():

    bin1 = bin_manager.Bin(contigs={"1", "2", "3"}, origin="A", name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "2", "4"}, origin="B", name="bin2")

    for b in [bin1, bin2]:
        b.completeness = 100
        b.contamination = 0
    
    G = nx.Graph()
    G.add_edge(bin1, bin2)

    return G


def test_get_intersection_bins(simple_bin_graph):

    intersec_bins = bin_manager.get_intersection_bins(simple_bin_graph)
    
    assert len(intersec_bins) == 1
    intersec_bin = intersec_bins.pop()

    assert intersec_bin.contigs == {"1", "2"}

def test_get_difference_bins(simple_bin_graph):
    
    difference_bins = bin_manager.get_difference_bins(simple_bin_graph)
    
    expected_bin1 = bin_manager.Bin(contigs={"3"}, origin="D", name="1")
    expected_bin2 = bin_manager.Bin(contigs={"4"}, origin="D", name="2")

    assert len(difference_bins) == 2
    assert difference_bins == {expected_bin1,expected_bin2}


def test_get_union_bins(simple_bin_graph):
    
    u_bins = bin_manager.get_union_bins(simple_bin_graph)
    
    expected_bin1 = bin_manager.Bin(contigs={"1", "2", "3", "4"}, origin="U", name="1")

    assert len(u_bins) == 1
    assert u_bins == {expected_bin1}




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
    result_bins = bin_manager.get_bins_from_contig2bin_table(str(test_table_path), set_name)

    # Validate the result
    assert len(result_bins) == 2  # Check if the correct number of bins are created

    # Define expected bins based on the test table content
    expected_bins = [
            bin_manager.Bin(contigs={"contig1", "contig2"}, origin="A", name="bin1"),
            bin_manager.Bin(contigs={"contig3"}, origin="A", name="bin2")
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
            "contig4\tbinA"
        ]
    }

    # Create temporary files for contig-to-bin tables
    for name, content in test_tables.items():
        table_path = tmp_path / f"test_{name}_contig2bin_table.txt"
        table_path.write_text("\n".join(content))

    # Call the function to parse contig-to-bin tables
    result_bin_dict = bin_manager.parse_contig2bin_tables({name: str(tmp_path / f"test_{name}_contig2bin_table.txt") for name in test_tables})

    # Validate the result
    assert len(result_bin_dict) == len(test_tables)  # Check if the number of bins matches the number of tables

    # Define expected Bin objects based on the test tables
    expected_bins = {
        "set1": [
            bin_manager.Bin(contigs={"contig1", "contig2"}, origin="set1", name="bin1"),
            bin_manager.Bin(contigs={"contig3"}, origin="set1", name="bin2"),
        ],
        "set2": [
            bin_manager.Bin(contigs={"contig3", "contig4"}, origin="set2", name="binA"),
        ]
    }

    # Compare expected bins with the result
    for name, expected in expected_bins.items():
        assert name in result_bin_dict
        assert len(result_bin_dict[name]) == len(expected)
        for result_bin in result_bin_dict[name]:
            assert result_bin in expected


def test_parse_contig2bin_tables_with_duplicated_bins(tmp_path, caplog):
    # Create temporary contig-to-bin tables for testing
    test_tables = {
        "set1": [
            "# Sample contig-to-bin table for bin1",
            "contig1\tbin1",
            "contig2\tbin1",
            "contig3\tbin2",
            "contig3\tbin3",
        ]
    }

    # Create temporary files for contig-to-bin tables
    for name, content in test_tables.items():
        table_path = tmp_path / f"test_{name}_contig2bin_table.txt"
        table_path.write_text("\n".join(content))

    # Call the function to parse contig-to-bin tables
    bin_manager.parse_contig2bin_tables({name: str(tmp_path / f"test_{name}_contig2bin_table.txt") for name in test_tables})
    expected_log_message = ('2 bins with identical contig compositions detected in bin set "set1". '
                           'These bins were merged to ensure uniqueness.')
    assert expected_log_message in caplog.text 


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

    bins = bin_manager.get_bins_from_directory(Path(bin_dir), set_name, fasta_extensions={'.fasta'})

    assert len(bins) == 2  # Ensure that the correct number of Bin objects is returned

    # Check if the Bin objects are created with the correct contigs, set name, and bin names
    assert isinstance(bins[0], bin_manager.Bin)
    assert isinstance(bins[1], bin_manager.Bin)
    assert bins[1].contigs in [{"contig1", "contig2"}, {"contig3", "contig4"}]
    assert bins[0].contigs in [{"contig1", "contig2"}, {"contig3", "contig4"}]
    assert bins[0].origin == {set_name}
    assert bins[1].origin == {set_name}
    assert bins[1].name in ["bin2.fasta", "bin1.fasta"]
    assert bins[0].name in ["bin2.fasta", "bin1.fasta"]

def test_get_bins_from_directory_no_files(tmpdir):
    bin_dir = Path(tmpdir.mkdir("empty_bins"))
    set_name = "EmptySet"

    bins = bin_manager.get_bins_from_directory(bin_dir, set_name, fasta_extensions={'.fasta'})

    assert len(bins) == 0  # Ensure that no Bin objects are returned for an empty directory

def test_get_bins_from_directory_no_wrong_extensions(create_temp_bin_files):
    bin_dir = Path(create_temp_bin_files)
    set_name = "TestSet"

    bins = bin_manager.get_bins_from_directory(bin_dir, set_name, fasta_extensions={'.fna'})

    assert len(bins) == 0  # Ensure that no Bin objects are returned for an empty directory





def test_parse_bin_directories(create_temp_bin_directories):
    set_name_to_bin_dir = create_temp_bin_directories

    bins = bin_manager.parse_bin_directories(set_name_to_bin_dir, fasta_extensions={'.fasta'})

    assert len(bins) == 2  # Ensure that the correct number of bin directories is parsed

    # Check if the Bin objects are created with the correct contigs, set name, and bin names
    assert isinstance(list(bins["set1"])[0], bin_manager.Bin)
    assert isinstance(list(bins["set2"])[0], bin_manager.Bin)

    assert len(bins["set2"]) == 1
    assert len(bins["set1"]) == 2


def test_get_contigs_in_bin_sets(example_bin_set1, example_bin_set2, caplog):
    """
    Test the get_contigs_in_bin_sets function for correct behavior.
    
    :param mock_bins: The mock_bins fixture providing test bin data.
    :param caplog: The pytest caplog fixture to capture logging output.
    """

    bin_set_name_to_bins = {"set1":example_bin_set1,
                            "set2":example_bin_set2}

    # Test the function with valid data
    with caplog.at_level(logging.WARNING):
        result = bin_manager.get_contigs_in_bin_sets(bin_set_name_to_bins)
    
    # Expected unique contigs
    expected_contigs = {"1", "2", "3", "4", "5"}
    
    # Check if the result matches expected contigs
    assert result == expected_contigs, "The returned set of contigs is incorrect."
    
def test_get_contigs_in_bin_sets_with_duplicated_warning(example_bin_set1, caplog):

    bin1 = bin_manager.Bin(contigs={"contig1", "2"}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={"contig1"}, origin="test1", name="binA")

    bin_set_name_to_bins = {
                            "set1":example_bin_set1,
                            "set_dup":{bin1, bin2},
                            }

    # Test the function with valid data
    with caplog.at_level(logging.WARNING):
        result = bin_manager.get_contigs_in_bin_sets(bin_set_name_to_bins)
    
    # Expected unique contigs
    expected_contigs = {"1", "2", "3", "4", "5", "contig1"}
    
    # Check if the result matches expected contigs
    assert result == expected_contigs, "The returned set of contigs is incorrect."

    # Check for expected warnings about duplicate contigs
    duplicate_warning = "Bin set 'set_dup' contains 1 duplicated contigs. Details: contig1 (found 2 times)"
    assert duplicate_warning in caplog.text, "The warning for duplicate contigs was not logged correctly."
