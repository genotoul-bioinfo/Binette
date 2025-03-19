from binette import cds
import pytest
import pyrodigal

from pathlib import Path
from unittest.mock import mock_open, patch

import gzip


class MockContig:
    def __init__(self, name, seq):
        self.seq = seq
        self.name = name


@pytest.fixture
def contig1():
    contig = MockContig(
        name="contig1",
        seq="ATGAGCATCAGGGGAGTAGGAGGATGCAACGGGAATAGTCGAATCCCTTCTCATAATGGGGATGGATCGAATCGCAGAAGTCAAAATACGAAGGGTAATAATAAAGTTGAAGATCGAGTTTGT",
    )
    return contig.name, contig.seq


@pytest.fixture
def contig2():
    contig = MockContig(
        name="contig2",
        seq="TTGGTCGTATGACTGATAATTTCTCAGACATTGAAAACTTTAATGAAATTTTCAACAGAAAACCTGCTTTACAATTTCGTTTTTA",
    )
    return contig.name, contig.seq


@pytest.fixture
def orf_finder():
    try:
        # for version >=3 of pyrodigal
        orf_finder = pyrodigal.GeneFinder(meta="meta")
    except AttributeError:
        orf_finder = pyrodigal.OrfFinder(meta="meta")

    return orf_finder


# Predict open reading frames with Pyrodigal using 1 thread.
def test_predict_orf_with_1_thread(contig1, contig2):

    contigs_iterator = [contig1, contig2]
    outfaa = "output.fasta"
    threads = 1

    result = cds.predict(contigs_iterator, outfaa, threads)

    assert isinstance(result, dict)
    assert len(result) == 2
    assert "contig1" in result
    assert "contig2" in result
    assert isinstance(result["contig1"], list)
    assert isinstance(result["contig2"], list)
    assert len(result["contig1"]) == 1
    assert len(result["contig2"]) == 1
    assert isinstance(result["contig1"][0], str)
    assert isinstance(result["contig2"][0], str)


def test_predict_orf_with_multiple_threads(contig1, contig2):

    contigs_iterator = [contig1, contig2]
    outfaa = "output.fasta"
    threads = 4

    result = cds.predict(contigs_iterator, outfaa, threads)

    assert isinstance(result, dict)
    assert len(result) == 2
    assert "contig1" in result
    assert "contig2" in result
    assert isinstance(result["contig1"], list)
    assert isinstance(result["contig2"], list)
    assert len(result["contig1"]) == 1
    assert len(result["contig2"]) == 1
    assert isinstance(result["contig1"][0], str)
    assert isinstance(result["contig2"][0], str)


def test_predict_genes(contig1, orf_finder):
    name, seq = contig1
    result = cds.predict_genes(orf_finder.find_genes, name, seq)

    assert isinstance(result, tuple)
    assert len(result) == 2
    assert isinstance(result[0], str)
    assert result[0] == name
    assert result[0] == "contig1"


# Extract the contig name from a CDS name.
def test_extract_contig_name_from_cds_name():
    cds_name = "contig1_gene1"

    result = cds.get_contig_from_cds_name(cds_name)

    assert isinstance(result, str)
    assert result == "contig1"


def test_write_faa(contig1, orf_finder):
    name, seq = contig1
    predicted_genes = orf_finder.find_genes(seq)
    contig_name = "contig"
    output_file = "tests/tmp_file.faa.gz"

    cds.write_faa(output_file, [(contig_name, predicted_genes)])

    # Check if the file was created and first seq starts
    # with the contig name as expected
    assert Path(output_file).exists()
    with gzip.open(output_file, "rt") as f:
        assert f.read().startswith(f">{contig_name}")


def test_parse_faa_file(tmp_path):
    # Mock a FASTA file of protein sequences
    # at least one protein sequence to not triger the error
    fasta_content = (
        ">contig1_gene1\n"
        "MPPPAOSKNSKSS\n"
        ">contig1_gene2\n"
        "CCCCCCCCCCC\n"
        ">contig2_gene1\n"
        "TTTTTTTTTTTT\n"
    )
    faa_file = tmp_path / "mock_file.faa"
    faa_file.write_text(fasta_content)

    # Call the function
    result = cds.parse_faa_file(faa_file)

    # Check if the output matches the expected dictionary
    expected_result = {
        "contig1": ["MPPPAOSKNSKSS", "CCCCCCCCCCC"],
        "contig2": ["TTTTTTTTTTTT"],
    }
    assert result == expected_result


def test_parse_faa_file_raises_error_for_dna(tmp_path):
    # Mock a DNA FASTA file
    fasta_content = (
        ">contig1_gene1\n"
        "AAAAAAAAAAA\n"
        ">contig1_gene2\n"
        "CCCCCCCCCCC\n"
        ">contig2_gene1\n"
        "TTTTTTTTTTTT\n"
    )
    fna_file = tmp_path / "mock_file.fna"
    fna_file.write_text(fasta_content)

    # Check that ValueError is raised when DNA sequences are encountered
    with pytest.raises(ValueError):
        cds.parse_faa_file(fna_file)


def test_get_aa_composition():

    genes = ["AAAA", "CCCC", "TTTT", "GGGG"]

    result = cds.get_aa_composition(genes)

    assert dict(result) == {"A": 4, "C": 4, "T": 4, "G": 4}


def test_get_contig_cds_metadata_flat():

    contig_to_genes = {"c1": ["AAAA", "GGGG", "CCCC"], "c2": ["TTTT", "CCCC"]}

    contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length = (
        cds.get_contig_cds_metadata_flat(contig_to_genes)
    )

    assert contig_to_cds_count == {"c1": 3, "c2": 2}
    assert contig_to_aa_counter == {
        "c1": {"A": 4, "G": 4, "C": 4},
        "c2": {"C": 4, "T": 4},
    }
    assert contig_to_aa_length == {"c1": 12, "c2": 8}


def test_get_contig_cds_metadata():

    contig_to_genes = {"c1": ["AAAA", "GGGG", "CCCC"], "c2": ["TTTT", "CCCC"]}

    contig_metadata = cds.get_contig_cds_metadata(contig_to_genes, 1)

    assert contig_metadata["contig_to_cds_count"] == {"c1": 3, "c2": 2}
    assert contig_metadata["contig_to_aa_counter"] == {
        "c1": {"A": 4, "G": 4, "C": 4},
        "c2": {"C": 4, "T": 4},
    }
    assert contig_metadata["contig_to_aa_length"] == {"c1": 12, "c2": 8}


# Test function
def test_is_nucleic_acid():
    # Valid DNA sequence
    assert cds.is_nucleic_acid("ATCG") is True
    assert cds.is_nucleic_acid("ATCNNNNNG") is True  # N can be found in DNA seq
    # Valid RNA sequence
    assert cds.is_nucleic_acid("AUGCAUGC") is True

    # Mixed case
    assert cds.is_nucleic_acid("AtCg") is True

    # Invalid sequence (contains characters not part of DNA or RNA)
    assert cds.is_nucleic_acid("ATCX") is False  # 'X' is not a valid base
    assert cds.is_nucleic_acid("AUG#C") is False  # '#' is not a valid base

    # Amino acid sequence
    assert cds.is_nucleic_acid("MSIRGVGGNGNSR") is False  # Numbers are invalid


def test_filter_faa_file_basic(tmp_path):
    """Test basic functionality of filtering a FASTA file."""
    # Create input data
    input_faa = tmp_path / "input.faa"
    filtered_faa = tmp_path / "filtered.faa"

    input_faa.write_text(
        ">contig1_gene1\nATGCGT\n" ">contig2_gene1\nATGCCG\n" ">contig3_gene1\nATGAAA\n"
    )

    # Contigs to keep
    contigs_to_keep = {"contig1", "contig3"}

    # Run the function
    cds.filter_faa_file(contigs_to_keep, input_faa, filtered_faa)

    # Check the filtered output
    expected_output = ">contig1_gene1\nATGCGT\n>contig3_gene1\nATGAAA\n"
    assert filtered_faa.read_text() == expected_output


def test_filter_faa_file_gz_output(tmp_path):
    """Test filtering with gzipped output."""
    input_faa = tmp_path / "input.faa"
    filtered_faa = tmp_path / "filtered.faa.gz"

    input_faa.write_text(
        ">contig1_gene1\nATGCGT\n" ">contig2_gene1\nATGCCG\n" ">contig3_gene1\nATGAAA\n"
    )

    # Contigs to keep
    contigs_to_keep = {"contig2"}

    # Run the function
    cds.filter_faa_file(contigs_to_keep, input_faa, filtered_faa)

    # Check the filtered output
    with gzip.open(filtered_faa, "rt") as f:
        filtered_content = f.read()
    expected_output = ">contig2_gene1\nATGCCG\n"
    assert filtered_content == expected_output


def test_filter_faa_file_no_matching_contigs(tmp_path):
    """Test filtering when no contigs match the input list."""
    input_faa = tmp_path / "input.faa"
    filtered_faa = tmp_path / "filtered_no_match.faa"

    input_faa.write_text(
        ">contig1_gene1\nATGCGT\n" ">contig2_gene1\nATGCCG\n" ">contig3_gene1\nATGAAA\n"
    )

    # Contigs to keep
    contigs_to_keep = {"contig4", "contig5"}

    # Run the function
    cds.filter_faa_file(contigs_to_keep, input_faa, filtered_faa)

    # Check the output file is empty
    assert filtered_faa.read_text() == ""
