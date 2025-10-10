import logging
import os
import sys
from collections import Counter
from pathlib import Path
from unittest.mock import patch

import pytest
from pyroaring import BitMap

from binette.bin_manager import Bin
from binette.main import (
    log_selected_bin_info,
    main,
    manage_protein_alignement,
    parse_input_files,
)
from tests.unit_tests.bin_manager_test import (  # noqa: F401
    create_temp_bin_directories,
    create_temp_bin_files,
)


@pytest.fixture
def test_environment(tmp_path: Path):
    """
    Fixture to set up a test environment with required directories and files.
    """
    folder1 = tmp_path / "folder1"
    folder2 = tmp_path / "folder2"
    contigs_file = tmp_path / "contigs.fasta"

    folder1.mkdir()
    folder2.mkdir()
    contigs_file.write_text(">contig1\nATCG")  # Sample content for the FASTA file

    return folder1, folder2, contigs_file


@pytest.fixture
def bins():
    b1 = Bin(contigs=BitMap({1}), origin="set1", name="bin1")
    b2 = Bin(contigs=BitMap({3}), origin="set1", name="bin2")
    b3 = Bin(contigs=BitMap({3, 2}), origin="set1", name="bin3")

    b1.add_quality(100, 0, 0)
    b2.add_quality(95, 10, 0)
    b3.add_quality(70, 20, 0)

    return [b1, b2, b3]


def test_log_selected_bin_info(caplog, bins):
    caplog.set_level(logging.INFO)

    hq_min_completeness = 85
    hq_max_conta = 15

    # Call the function
    log_selected_bin_info(bins, hq_min_completeness, hq_max_conta)

    # Check if the logs contain expected messages
    expected_logs = "2/3 selected bins have high quality (completeness >= 85 and contamination <= 15)"

    assert expected_logs in caplog.text


def test_manage_protein_alignement_resume(tmp_path):
    # Create temporary directories and files for testing

    faa_file = tmp_path / "proteins.faa"
    faa_file_content = (
        ">contig1_1\nMCGT\n>contig2_1\nTGCA\n>contig2_2\nAAAA\n>contig3_1\nCCCC\n"
    )

    faa_file.write_text(faa_file_content)

    contig_to_kegg_id = {
        "contig1": Counter({"K12345": 1, "K67890": 1}),
        "contig2": Counter({"K23456": 1}),
    }

    with patch("binette.diamond.get_contig_to_kegg_id", return_value=contig_to_kegg_id):
        # Call the function

        # Run the function with test data
        contig_to_kegg_counter, contig_to_genes, _ = manage_protein_alignement(
            faa_file=Path(faa_file),
            contigs_fasta=Path("contigs_fasta"),
            contigs_in_bins=set(("contig1", "contig2", "contig3")),
            diamond_result_file=Path("diamond_result_file"),
            checkm2_db=None,
            threads=1,
            use_existing_protein_file=True,
            resume_diamond=True,
            low_mem=False,
        )

    # Assertions to check the function output or file existence
    assert isinstance(contig_to_genes, dict)
    assert isinstance(contig_to_kegg_counter, dict)
    assert len(contig_to_genes) == 3


def test_manage_protein_alignement_not_resume(tmpdir, tmp_path):
    # Create temporary directories and files for testing

    faa_file = tmp_path / "proteins.faa"
    faa_file_content = ">contig1_1\nMLKPACGT\n>contig2_1\nMMMKPTGCA\n>contig2_2\nMMMAAAA\n>contig3_1\nMLPALP\n"

    faa_file.write_text(faa_file_content)

    contigs_fasta = os.path.join(str(tmpdir), "contigs.fasta")
    diamond_result_file = os.path.join(str(tmpdir), "diamond_results.tsv")

    contig_to_kegg_id = {
        "contig1": Counter({"K12345": 1, "K67890": 1}),
        "contig2": Counter({"K23456": 1}),
    }

    with (
        patch("binette.diamond.get_contig_to_kegg_id", return_value=contig_to_kegg_id),
        patch("binette.diamond.run", return_value=None),
    ):
        # Call the function

        contig_to_kegg_counter, contig_to_genes, _ = manage_protein_alignement(
            faa_file=Path(faa_file),
            contigs_fasta=Path(contigs_fasta),
            contigs_in_bins=set(("contig1", "contig2", "contig3")),
            diamond_result_file=Path(diamond_result_file),
            checkm2_db=None,
            threads=1,
            use_existing_protein_file=True,
            resume_diamond=True,
            low_mem=False,
        )

    # Assertions to check the function output or file existence
    assert isinstance(contig_to_genes, dict)
    assert isinstance(contig_to_kegg_counter, dict)
    assert len(contig_to_genes) == 3


def test_parse_input_files_with_contig2bin_tables(tmp_path):
    bin_set1 = tmp_path / "bin_set1.tsv"
    bin_set1.write_text("contig1\tbin1A\ncontig2\tbin1B\n")
    bin_set2 = tmp_path / "bin_set2.tsv"
    bin_set2.write_text("contig3\tbin2A\ncontig4\ttbin2B\n")

    fasta_file = tmp_path / "assembly.fasta"
    fasta_file_content = ">contig1\nACGT\n>contig2\nTGCA\n>contig3\nAAAA\n>contig4\nCCCC\n>contig5\nCGTCGCT\n"
    fasta_file.write_text(fasta_file_content)

    # Call the function and capture the return values
    (
        contig_key_to_bin,
        contigs_in_bins,
        contig_id_to_length,
        contig_to_index,
    ) = parse_input_files(None, [bin_set1, bin_set2], fasta_file, tmp_path)

    # # Perform assertions on the returned values
    assert isinstance(contig_key_to_bin, dict)
    assert isinstance(contigs_in_bins, list)
    assert isinstance(contig_id_to_length, dict)

    assert len(contig_key_to_bin) == 4
    assert set(contigs_in_bins) == {"contig1", "contig2", "contig3", "contig4"}
    assert len(contig_id_to_length) == 4


def test_parse_input_files_with_contig2bin_tables_with_unknown_contig(tmp_path):
    bin_set3 = tmp_path / "bin_set3.tsv"
    bin_set3.write_text("contig3\tbin3A\ncontig44\ttbin3B\n")
    fasta_file = tmp_path / "assembly.fasta"
    fasta_file_content = ">contig1\nACGT\n>contig2\nTGCA\n>contig3\nAAAA\n>contig4\nCCCC\n>contig5\nCGTCGCT\n"
    fasta_file.write_text(fasta_file_content)

    with pytest.raises(ValueError):
        parse_input_files(None, [bin_set3], fasta_file, tmp_path)


def test_parse_input_files_bin_dirs(create_temp_bin_directories, tmp_path):
    bin_dirs = [Path(d) for d in create_temp_bin_directories.values()]

    contig2bin_tables = []

    # Create temporary directories and files for testing

    fasta_file = tmp_path / "assembly.fasta"
    fasta_file_content = ">contig1\nACGT\n>contig2\nTGCA\n>contig3\nAAAA\n>contig4\nCCCC\n>contig5\nCGTCGCT\n"
    fasta_file.write_text(fasta_file_content)

    # Call the function and capture the return values
    (
        contig_key_to_bin,
        contigs_in_bins,
        contig_id_to_length,
        contig_to_index,
    ) = parse_input_files(bin_dirs, contig2bin_tables, fasta_file)

    # # Perform assertions on the returned values
    assert isinstance(contig_key_to_bin, dict)
    assert isinstance(contigs_in_bins, list)
    assert isinstance(contig_id_to_length, dict)

    assert len(contig_key_to_bin) == 3
    assert set(contigs_in_bins) == {
        "contig1",
        "contig2",
        "contig3",
        "contig4",
        "contig5",
    }
    assert len(contig_id_to_length) == 5


# @patch('diamond.run')
def test_manage_protein_alignment_no_resume(tmp_path):
    # Set up the input parameters
    faa_file = Path("test.faa")
    contigs_fasta = Path("test.fasta")
    contigs_in_bins = {"bin1": ["contig1"]}
    diamond_result_file = Path("test_diamond_result.txt")
    checkm2_db = tmp_path / "checkm2_db"
    with open(checkm2_db, "w"):
        pass
    threads = 4
    resume = False
    low_mem = False

    # Mock the necessary functions
    with (
        patch("pyfastx.Fastx") as mock_pyfastx_Fastx,
        patch("binette.cds.predict") as mock_predict,
        patch("binette.diamond.get_checkm2_db"),
        patch("binette.diamond.run") as mock_diamond_run,
        patch(
            "binette.diamond.get_contig_to_kegg_id"
        ) as mock_diamond_get_contig_to_kegg_id,
    ):
        # Set the return value of the mocked functions
        mock_pyfastx_Fastx.return_value = [("contig1", "ATCG")]
        mock_predict.return_value = (
            {
                "contig1": ["gene1"],
            },
            {"contig1": 50},
        )

        # Call the function
        contig_to_kegg_counter, contig_to_genes, _ = manage_protein_alignement(
            faa_file,
            contigs_fasta,
            contigs_in_bins,
            diamond_result_file,
            checkm2_db,
            threads,
            resume,
            resume,
            low_mem,
        )

        # Assertions to check if functions were called
        mock_pyfastx_Fastx.assert_called_once()
        mock_predict.assert_called_once()
        mock_diamond_get_contig_to_kegg_id.assert_called_once()
        mock_diamond_run.assert_called_once_with(
            faa_file.as_posix(),
            diamond_result_file.as_posix(),
            checkm2_db.as_posix(),
            f"{os.path.splitext(diamond_result_file.as_posix())[0]}.log",
            threads,
            low_mem=low_mem,
        )


def test_main_resume_when_not_possible(monkeypatch, test_environment, tmp_path):
    # Define or mock the necessary inputs/arguments
    folder1, folder2, contigs_file = test_environment

    outdir = tmp_path / "results"
    # Mock sys.argv to use test_args
    test_args = [
        "-d",
        str(folder1),
        str(folder2),
        "-c",
        str(contigs_file),
        # ... more arguments as required ...
        "--debug",
        "--resume",
        "-o",
        outdir.as_posix(),
    ]
    monkeypatch.setattr(sys, "argv", ["your_script.py"] + test_args)

    # Call the main function
    with pytest.raises(FileNotFoundError):
        main()
