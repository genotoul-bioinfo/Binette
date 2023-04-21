"""
Unit tests for binette.

"""

import pytest

# from bin_manager import get_all_possible_combinations
# from . import bin_manager
from binette import bin_manager


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


def test_two_bin_intersection():
    bin1 = bin_manager.Bin(contigs={"1", "2", "e", "987"}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={"1", "e", "2", "33"}, origin="test2", name="binA")

    bin_intersection = bin1 & bin2

    assert bin_intersection == bin_manager.Bin({"1", "2", "e"}, "", "")


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

    union_bin = bin_manager.Bin(contigs={"13", "21", "1", "e", "2", "33"}, origin="", name="")

    assert bin1.union(bin2) == union_bin
    assert bin1.union(bin2).name == "bin1 | binA"


def test_bin_difference():
    bin1 = bin_manager.Bin(contigs={5, 1, 6, 7, 8}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={3, 6, 7}, origin="test2", name="bin2")
    bin3 = bin_manager.Bin(contigs={1, 2, 3, 4}, origin="test2", name="bin3")

    diff_bin1_23 = bin_manager.Bin(contigs={5, 8}, origin="", name="")
    diff_bin1_2 = bin_manager.Bin(contigs={5, 1, 8}, origin="", name="")

    assert bin1.difference(bin2, bin3) == diff_bin1_23
    assert bin1.difference(bin2) == diff_bin1_2
    assert bin1.difference(bin2, bin3).name == "bin1 - bin2 - bin3"


def test_bin_intersection():
    bin1 = bin_manager.Bin(contigs={5, 1, 6, 7, 8}, origin="test1", name="bin1")
    bin2 = bin_manager.Bin(contigs={3, 6, 7}, origin="test2", name="bin2")
    bin3 = bin_manager.Bin(contigs={1, 2, 3, 4, 7}, origin="test2", name="bin3")

    inter_bin123 = bin_manager.Bin(contigs={7}, origin="", name="")
    iner_bin1_2 = bin_manager.Bin(contigs={7, 6}, origin="", name="")

    assert bin1.intersection(bin2, bin3) == inter_bin123
    assert bin1.intersection(bin2) == iner_bin1_2
    assert bin1.intersection(bin2, bin3).name == "bin1 & bin2 & bin3"


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
