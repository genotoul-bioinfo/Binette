'''
Unit tests for bin_refinement.

'''

import pytest 

# from bin_manager import get_all_possible_combinations
import bin_manager

def test_get_all_possible_combinations():
    input_list = ['2', "3", "4"]
    expected_list = [("2", "3"), 
                        ("2", "4"), 
                        ("3", "4"),
                        ("2", "3", "4")]

    assert list(bin_manager.get_all_possible_combinations(input_list)) == expected_list


@pytest.fixture
def example_bin_set1():
    bin1 = bin_manager.Bin(contigs={"1", "2"}, origin='test1', name="bin1")
    bin2 = bin_manager.Bin(contigs={"3","4"}, origin='test1', name="bin2")
    bin3 = bin_manager.Bin(contigs={"5"}, origin='test1', name="bin2")
    return {bin1, bin2, bin3}


def test_bin_eq_true():
    bin1 = bin_manager.Bin(contigs={"1", "2", "e"}, origin='test1', name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2"}, origin='test2', name="binA")

    assert bin1 == bin2

def test_bin_eq_false():
    bin1 = bin_manager.Bin(contigs={"1", "2", "e"}, origin='test1', name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin='test2', name="binA")

    assert bin1 != bin2

def test_in_for_bin_list():
    bin1 = bin_manager.Bin(contigs={"1", "2", "e"}, origin='test1', name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin='test2', name="binA")
    bin3 = bin_manager.Bin(contigs={"4", "R"}, origin='test2', name="binA")

    bins = [bin1, bin2]

    assert bin1 in bins 
    assert bin2 in bins
    assert bin3 not in bins

def test_two_bin_intersection():
    bin1 = bin_manager.Bin(contigs={"1", "2", "e", "987"}, origin='test1', name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin='test2', name="binA")
    
    bin_intersection = bin1 & bin2

    assert bin_intersection == bin_manager.Bin({"1", "2", "e"}, "", "")

def test_multiple_bins_intersection():
    bin1 = bin_manager.Bin(contigs={"1", "2", "e", "987"}, origin='test1', name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin='test2', name="binA")
    bin3 = bin_manager.Bin(contigs={"1", "123", "2", "33"}, origin='test2', name="binA")
    
    bin_intersection = bin1.intersection(bin2, bin3)

    assert bin_intersection == bin_manager.Bin({"1", "2"}, "", "")



def test_bin_set_has_no_duplicate(example_bin_set1):
    set1_twice = list(example_bin_set1) + list(example_bin_set1) 

    assert example_bin_set1 == set(set1_twice)

def test_bin_overlap_true():
    bin1 = bin_manager.Bin(contigs={"1", "2", "e"}, origin='test1', name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin='test2', name="binA")

    assert bin1.overlaps_with(bin2)

    assert bin2.overlaps_with(bin1)



def test_bin_overlap_flase():
    bin1 = bin_manager.Bin(contigs={"13", "21", "ef"}, origin='test1', name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin='test2', name="binA")

    assert not bin1.overlaps_with(bin2)
    
    