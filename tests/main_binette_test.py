import pytest
import logging
from binette.main import (
    log_selected_bin_info,
    manage_protein_alignement,
    parse_input_files,
    parse_arguments,
    init_logging,
    main,
    UniqueStore,
    is_valid_file,
)
from binette.bin_manager import Bin
from binette import diamond, contig_manager, cds
import os
import sys
from unittest.mock import patch, MagicMock

from collections import Counter
from tests.bin_manager_test import create_temp_bin_directories, create_temp_bin_files
from argparse import ArgumentParser
from pathlib import Path
from pyroaring import BitMap


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
    expected_logs = "2/3 selected bins have a high quality (completeness >= 85 and contamination <= 15)."

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
        contig_to_kegg_counter, contig_to_genes = manage_protein_alignement(
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

    contig_to_length = {"contig1": 40, "contig2": 80, "contig3": 20}

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

        contig_to_kegg_counter, contig_to_genes = manage_protein_alignement(
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


def test_argument_used_once():
    # Test UniqueStore class
    parser = ArgumentParser(description="Test parser")
    parser.add_argument("--example", action=UniqueStore, help="Example argument")
    args = parser.parse_args(["--example", "value"])
    assert args.example == "value"


def test_argument_used_multiple_times():
    # Test UniqueStore class
    parser = ArgumentParser(description="Test parser")
    parser.add_argument("--example", action=UniqueStore, help="Example argument")
    with pytest.raises(SystemExit):
        parser.parse_args(["--example", "value", "--example", "value2"])


def test_parse_arguments_required_arguments(test_environment):
    """
    Test parsing when only required arguments are provided.
    Ensure that input arguments exist before parsing.
    """
    # Create temporary directories and files
    folder1, folder2, contigs_file = test_environment

    # Parse arguments with existing files and directories
    args = parse_arguments(["-d", str(folder1), str(folder2), "-c", str(contigs_file)])

    # Assert that the parsed arguments match the expected paths
    assert args.bin_dirs == [folder1, folder2]
    assert args.contigs == contigs_file


def test_parse_arguments_optional_arguments(test_environment):
    # Test when required and optional arguments are provided

    # Create temporary directories and files
    folder1, folder2, contigs_file = test_environment

    # Parse arguments with existing files and directories
    args = parse_arguments(
        [
            "-d",
            str(folder1),
            str(folder2),
            "-c",
            str(contigs_file),
            "--threads",
            "4",
            "--outdir",
            "output",
        ]
    )
    assert args.bin_dirs == [folder1, folder2]
    assert args.contigs == contigs_file
    assert args.threads == 4
    assert args.outdir == Path("output")


def test_parse_arguments_invalid_arguments():
    # Test when invalid arguments are provided
    with pytest.raises(SystemExit):
        # In this case, required arguments are missing
        parse_arguments(["-t", "4"])


def test_parse_arguments_help():
    # Test the help message
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parse_arguments(["-h"])
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


def test_init_logging_command_line(caplog):

    caplog.set_level(logging.INFO)

    init_logging(verbose=True, debug=False)
    expected_log_message = f'command line: {" ".join(sys.argv)}'
    # Check if the log message is present in the log records

    assert expected_log_message in caplog.text


# @patch('diamond.run')
def test_manage_protein_alignment_no_resume(tmp_path):
    # Set up the input parameters
    faa_file = Path("test.faa")
    contigs_fasta = Path("test.fasta")
    contig_to_length = {"contig1": [1000]}
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
        patch("binette.diamond.get_checkm2_db") as mock_get_checkm2_db,
        patch("binette.diamond.run") as mock_diamond_run,
        patch(
            "binette.diamond.get_contig_to_kegg_id"
        ) as mock_diamond_get_contig_to_kegg_id,
    ):

        # Set the return value of the mocked functions
        mock_pyfastx_Fastx.return_value = [("contig1", "ATCG")]
        mock_predict.return_value = {"contig1": ["gene1"]}

        # Call the function
        contig_to_kegg_counter, contig_to_genes = manage_protein_alignement(
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


def test_main_resume_when_not_possible(monkeypatch, test_environment):
    # Define or mock the necessary inputs/arguments
    folder1, folder2, contigs_file = test_environment

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
    ]
    monkeypatch.setattr(sys, "argv", ["your_script.py"] + test_args)

    # Call the main function
    with pytest.raises(FileNotFoundError):
        main()


def test_is_valid_file_existing_file(tmp_path: Path):
    """Test is_valid_file with a file that exists."""
    # Create a temporary file
    test_file = tmp_path / "test_file.txt"
    test_file.write_text("Sample content")

    parser = ArgumentParser()

    # Assert that the function correctly returns the file path
    result = is_valid_file(parser, str(test_file))
    assert result == test_file


def test_is_valid_file_non_existing_file():
    """Test is_valid_file with a file that does not exist."""
    parser = ArgumentParser()
    non_existing_file = "non_existing_file.txt"

    # Expect the function to call parser.error, which will raise a SystemExit exception
    with pytest.raises(SystemExit):
        is_valid_file(parser, non_existing_file)
