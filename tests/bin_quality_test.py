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
     get_bins_metadata_df)

from checkm2 import keggData, modelPostprocessing, modelProcessing

from unittest.mock import Mock, patch, MagicMock

def test_compute_N50():
    assert bin_quality.compute_N50([50]) == 50
    assert bin_quality.compute_N50([0]) == 0
    assert bin_quality.compute_N50([]) == 0
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
    


class Bin:
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
    bins = [
        Bin(1, ['contig1', 'contig3']),
        Bin(2, ['contig2'])
    ]

    contig_to_cds_count = {'contig1': 10, 'contig2': 45, 'contig3': 20, 'contig4': 25}
    contig_to_aa_counter = {'contig1': Counter({'A': 5, 'D': 10}), 'contig2': Counter({'G': 8, 'V': 12, 'T': 2}), 
                            'contig3': Counter({'D': 8, 'Y': 12})}
    contig_to_aa_length = {'contig1': 1000, 'contig2': 1500, 'contig3': 2000, 'contig4': 2500}

    # Call the function
    result_df = bin_quality.get_bins_metadata_df(bins, contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length)

    # Define expected values based on the provided input
    expected_columns = [
        'Name', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
        'AALength', 'CDS'
    ]
    expected_index = [1, 2]

    expected_values = [
        [1, 5, 0, 18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 3000, 30],
        [2, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 12, 0, 0, 1500, 45]
    ]

    # Check if the generated DataFrame matches the expected DataFrame
    assert result_df.columns.tolist() == expected_columns
    assert result_df.index.tolist() == expected_index
    assert result_df.values.tolist() == expected_values




def test_get_diamond_feature_per_bin_df():
    # Mock input data
    bins = [
        Bin(1, ['contig1', 'contig2']),
        Bin(2, ['contig3', 'contig4'])
    ]

    contig_to_kegg_counter = {
        'contig1': Counter({'K01810': 5, 'K15916': 7}),
        'contig2': Counter({'K01810': 10}),
        'contig3': Counter({'K00918': 8}),
    }

    # Call the function
    result_df, default_ko_count = bin_quality.get_diamond_feature_per_bin_df(bins, contig_to_kegg_counter)

    expected_index = [1, 2]
    
    
    assert result_df.index.tolist() == expected_index
    assert result_df.loc[1,"K01810"] == 15 # in bin1 from contig 1 and 2
    assert result_df.loc[1,"K15916"] == 7 # in bin1 from contig 1 
    assert result_df.loc[2,"K01810"] == 0 # this ko is not in any contig of bin 2
    assert result_df.loc[2,"K00918"] == 8 # in bin2 from contig 3



def test_add_bin_size_and_N50():
    # Mock input data
    bins = [
        Bin(1, ['contig1', 'contig2']),
        Bin(2, ['contig3'])
    ]

    contig_to_size = {
        'contig1': 1000,
        'contig2': 1500,
        'contig3': 2000,
    }

    # Call the function
    bin_quality.add_bin_size_and_N50(bins, contig_to_size)

    # Assertions to verify if add_length and add_N50 were called with the correct values
    assert bins[0].length == 2500
    assert bins[0].N50 == 1500 
    assert bins[1].length == 2000
    assert bins[1].N50 == 2000 




def mock_modelProcessor(thread):
    return "mock_modelProcessor"

def test_add_bin_metrics(monkeypatch):
    # Mock input data
    bins = [
        Bin(1, ['contig1', 'contig2']),
        Bin(2, ['contig3'])
    ]

    contig_info = {
        # Add mocked contig information here as needed
        "contig_to_kegg_counter":{},
        "contig_to_cds_count":{},
        "contig_to_aa_counter":{},
        "contig_to_aa_length":{},
        "contig_to_length":{},
    }

    contamination_weight = 0.5
    threads = 1
    

    monkeypatch.setattr(modelPostprocessing, "modelProcessor", mock_modelProcessor)


    # Mock the functions called within add_bin_metrics
    with patch('binette.bin_quality.add_bin_size_and_N50') as mock_add_bin_size_and_N50, \
         patch('binette.bin_quality.assess_bins_quality_by_chunk') as mock_assess_bins_quality_by_chunk:

        add_bin_metrics(bins, contig_info, contamination_weight, threads)

        # Assertions to check if functions were called with the expected arguments
        mock_add_bin_size_and_N50.assert_called_once_with(bins, contig_info["contig_to_length"])
        mock_assess_bins_quality_by_chunk.assert_called_once_with(
            bins,
            contig_info["contig_to_kegg_counter"],
            contig_info["contig_to_cds_count"],
            contig_info["contig_to_aa_counter"],
            contig_info["contig_to_aa_length"],
            contamination_weight,
            "mock_modelProcessor",  # Mocked postProcessor object
        )

def test_assess_bins_quality_by_chunk(monkeypatch):
    # Prepare input data for testing
    bins = [
        Bin(1, ['contig1', 'contig2']),
        Bin(2, ['contig3', 'contig4']),
        Bin(3, ['contig3', 'contig4'])
    ]

    contig_to_kegg_counter = {}
    contig_to_cds_count = {}
    contig_to_aa_counter = {}
    contig_to_aa_length = {}  
    contamination_weight = 0.5 

    # Mocking postProcessor object

    monkeypatch.setattr(modelPostprocessing, "modelProcessor", mock_modelProcessor)


    # Mock the functions called within add_bin_metrics
    with patch('binette.bin_quality.assess_bins_quality') as mock_assess_bins_quality:

        assess_bins_quality_by_chunk(
                                    bins,
                                    contig_to_kegg_counter,
                                    contig_to_cds_count,
                                    contig_to_aa_counter,
                                    contig_to_aa_length,
                                    contamination_weight,
                                    postProcessor=None,
                                    threads=1,
                                    chunk_size = 3 
                                )

        # Chunk size > number of bin so only one chunk
        mock_assess_bins_quality.assert_called_once_with(
            bins= set(bins),
            contig_to_kegg_counter= contig_to_kegg_counter,
            contig_to_cds_count=contig_to_cds_count,
            contig_to_aa_counter=contig_to_aa_counter,
            contig_to_aa_length=contig_to_aa_length,
            contamination_weight=contamination_weight,
            postProcessor=None,
            threads=1
        )

    # Mock the functions called within add_bin_metrics
    with patch('binette.bin_quality.assess_bins_quality') as mock_assess_bins_quality:

        assess_bins_quality_by_chunk(
                                    bins,
                                    contig_to_kegg_counter,
                                    contig_to_cds_count,
                                    contig_to_aa_counter,
                                    contig_to_aa_length,
                                    contamination_weight,
                                    postProcessor=None,
                                    threads=1,
                                    chunk_size = 2
                                )

        # Chunk size < number of bin so  2 chunks with [bin1,bin2] and [bin3]
        assert mock_assess_bins_quality.call_count == 2



from unittest.mock import patch, MagicMock
import numpy as np
import pandas as pd
from checkm2 import keggData, modelPostprocessing, modelProcessing


def test_assess_bins_quality():
    # Prepare mock input data for testing
    bins = [
        Bin(1, ['contig1', 'contig2']),
        Bin(2, ['contig3', 'contig4'])
    ]

    contig_to_kegg_counter = {}
    contig_to_cds_count = {}
    contig_to_aa_length = {}
    contig_to_aa_counter = {}
    contamination_weight = 0.5


    # Call the function being tested
    assess_bins_quality(
        bins,
        contig_to_kegg_counter,
        contig_to_cds_count,
        contig_to_aa_counter,
        contig_to_aa_length,
        contamination_weight
    )



    # Verify the expected calls to add_quality for each bin object
    for bin_obj in bins:
        assert bin_obj.completeness is not None
        assert bin_obj.contamination is not None
        assert bin_obj.score is not None
        assert bin_obj.score == bin_obj.completeness - bin_obj.contamination * contamination_weight