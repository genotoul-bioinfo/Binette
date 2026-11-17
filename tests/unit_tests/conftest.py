"""Fixtures shared by several unit test modules."""

from pathlib import Path

import pytest


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
