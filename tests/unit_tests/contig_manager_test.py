import pyfastx
import pytest

from binette import contig_manager


# Parses a valid FASTA file and returns a pyfastx.Fasta object.
def test_valid_fasta_file(tmp_path):
    # Arrange
    fasta_file = "tests/unit_tests/contigs.fasta"

    index_file = tmp_path / "contigs.fasta.fxi"
    result = contig_manager.parse_fasta_file(fasta_file, str(index_file))

    # Assert
    assert isinstance(result, pyfastx.Fasta)


# Parses an invalid FASTA file and raises an exception.
def test_invalid_fasta_file():
    #
    fasta_file = "tests/unit_tests/contig_manager_test.py"

    # Act and Assert
    with pytest.raises(RuntimeError):
        contig_manager.parse_fasta_file(fasta_file, "./index.fxi")


# The function returns a tuple containing two dictionaries.
def test_returns_tuple():
    contigs = ["contig1", "contig2", "contig3"]
    contig_to_index = contig_manager.make_contig_index(contigs)
    assert isinstance(contig_to_index, dict)
    assert contig_to_index == {"contig1": 0, "contig2": 1, "contig3": 2}


# The function returns a dictionary with the same number of items as the input dictionary.
def test_same_number_of_items():
    contig_to_index = {"contig1": 0, "contig2": 1, "contig3": 2}
    contig_to_info = {"contig1": "info1", "contig2": "info2", "contig3": "info3"}
    expected_result = {0: "info1", 1: "info2", 2: "info3"}

    result = contig_manager.apply_contig_index(contig_to_index, contig_to_info)

    assert len(result) == len(expected_result)
