from binette import cds
import pytest
import pyrodigal

from pathlib import Path
from unittest.mock import mock_open, patch

class MockContig:
    def __init__(self, name, seq):
        self.seq = seq
        self.name = name

@pytest.fixture
def contig1():
    contig =  MockContig(name="contig1", 
                          seq="ATGAGCATCAGGGGAGTAGGAGGATGCAACGGGAATAGTCGAATCCCTTCTCATAATGGGGATGGATCGAATCGCAGAAGTCAAAATACGAAGGGTAATAATAAAGTTGAAGATCGAGTTTGT")
    return contig

@pytest.fixture
def contig2():
    contig =  MockContig(name="contig2", 
                          seq="TTGGTCGTATGACTGATAATTTCTCAGACATTGAAAACTTTAATGAAATTTTCAACAGAAAACCTGCTTTACAATTTCGTTTTTA")
    return contig

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


    result = cds.predict_genes(orf_finder.find_genes, contig1)

    assert isinstance(result, tuple)
    assert len(result) == 2
    assert isinstance(result[0], str)
    assert result[0] == contig1.name
    assert result[0] == "contig1"


# Extract the contig name from a CDS name.
def test_extract_contig_name_from_cds_name():
    cds_name = "contig1_gene1"

    result = cds.get_contig_from_cds_name(cds_name)

    assert isinstance(result, str)
    assert result == "contig1"

def test_extract_contig_name_from_cds_name():
    cds_name = "contig1_gene1"

    result = cds.get_contig_from_cds_name(cds_name)

    assert isinstance(result, str)
    assert result == "contig1"


# Import the functions write_faa and parse_faa_file here


def test_write_faa(contig1, orf_finder):
    
    predicted_genes = orf_finder.find_genes(contig1.seq)
    contig_name = 'contig'
    output_file = "tests/tmp_file.faa"

    cds.write_faa(output_file, [(contig_name, predicted_genes)])

    # Check if the file was created and first seq starts
    # with the contig name as expected
    assert Path(output_file).exists()
    with open(output_file, "r") as f:
        assert f.read().startswith(f">{contig_name}")


def test_parse_faa_file(tmp_path):
    # Mock a FASTA file
    fasta_content = (
        ">contig1_gene1\n"
        "AAAAAAAAAAA\n"
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
        'contig1': ['AAAAAAAAAAA', 'CCCCCCCCCCC'],
        'contig2': ['TTTTTTTTTTTT']
    }
    assert result == expected_result

def test_get_aa_composition():

    genes = ['AAAA',
            "CCCC",
            "TTTT",
            "GGGG"]

    result = cds.get_aa_composition(genes)

    assert dict(result) == {'A': 4, 'C': 4, 'T': 4, 'G': 4}

def test_get_contig_cds_metadata_flat():

    contig_to_genes = {"c1":["AAAA", "GGGG", "CCCC"],
                       "c2":["TTTT", "CCCC"]}

    contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length = cds.get_contig_cds_metadata_flat(contig_to_genes)
    
    assert contig_to_cds_count == {"c1":3, "c2":2}
    assert contig_to_aa_counter == {"c1": {'A': 4, 'G': 4, "C":4} , "c2":{'C': 4, 'T': 4}}
    assert contig_to_aa_length == {"c1":12, "c2":8}

def test_get_contig_cds_metadata():

    contig_to_genes = {"c1":["AAAA", "GGGG", "CCCC"],
                       "c2":["TTTT", "CCCC"]}

    contig_metadata = cds.get_contig_cds_metadata(contig_to_genes, 1)
    
    assert contig_metadata['contig_to_cds_count'] == {"c1":3, "c2":2}
    assert contig_metadata['contig_to_aa_counter'] == {"c1": {'A': 4, 'G': 4, "C":4} , "c2":{'C': 4, 'T': 4}}
    assert contig_metadata['contig_to_aa_length'] == {"c1":12, "c2":8}