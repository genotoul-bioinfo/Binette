
import pytest
import logging
from binette.main import log_selected_bin_info, select_bins_and_write_them, manage_protein_alignement, parse_input_files, parse_arguments, init_logging, main, UniqueStore
from binette.bin_manager import Bin
from binette import diamond, contig_manager, cds
import os
import sys
from unittest.mock import patch, MagicMock

from collections import Counter
from tests.bin_manager_test import create_temp_bin_directories, create_temp_bin_files
from argparse import ArgumentParser
from pathlib import Path

@pytest.fixture
def bins():
    b1 = Bin(contigs={"contig1"}, origin="set1", name="bin1")
    b2 = Bin(contigs={"contig3"}, origin="set1", name="bin2")
    b3 = Bin(contigs={"contig3", "contig2"}, origin="set1", name="bin3")

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
    expected_logs ="2/3 selected bins have a high quality (completeness >= 85 and contamination <= 15)."


    assert expected_logs in caplog.text


def test_select_bins_and_write_them(tmp_path, tmpdir, bins):
    # Create temporary directories and files for testing
    outdir = tmpdir.mkdir("test_outdir")
    contigs_fasta = os.path.join(str(outdir), "contigs.fasta")
    final_bin_report = os.path.join(str(outdir), "final_bin_report.tsv")

    index_to_contig={"contig1":"contig1", "contig2": "contig2", "contig3":"contig3"}

    contigs_fasta = tmp_path / "contigs.fasta"
    contigs_fasta_content = (
    ">contig1\nACGT\n>contig2\nTGCA\n>contig3\nAAAA\n>contig4\nCCCC\n"
    )
    contigs_fasta.write_text(contigs_fasta_content)

    b1, b2, b3 = bins


    # Run the function with test data
    selected_bins = select_bins_and_write_them(
        set(bins), contigs_fasta, Path(final_bin_report), min_completeness=60, index_to_contig=index_to_contig, outdir=outdir, debug=True
    )

    # Assertions to check the function output or file existence
    assert isinstance(selected_bins, list)
    assert os.path.isfile(final_bin_report)
    assert selected_bins == bins[:2] # The third bin is overlapping with the second one and has a worse score so it is not selected.

    with open(outdir / f"final_bins/bin_{b1.id}.fa", "r") as bin1_file:
        assert bin1_file.read() == ">contig1\nACGT\n"

    with open(outdir / f"final_bins/bin_{b2.id}.fa", "r") as bin2_file:
        assert bin2_file.read() == ">contig3\nAAAA\n"

    assert not os.path.isfile(outdir / f"final_bins/bin_{b3.id}.fa")



def test_manage_protein_alignement_resume(tmp_path):
    # Create temporary directories and files for testing

    faa_file = tmp_path / "proteins.faa"
    faa_file_content = (
    ">contig1_1\nACGT\n>contig2_1\nTGCA\n>contig2_2\nAAAA\n>contig3_1\nCCCC\n"
    )

    contig_to_length={"contig1":40, "contig2":80, "contig3":20}

    faa_file.write_text(faa_file_content)

    contig_to_kegg_id = {
        "contig1": Counter({"K12345": 1, "K67890": 1}),
        "contig2": Counter({"K23456": 1})
    }


    with patch("binette.diamond.get_contig_to_kegg_id", return_value=contig_to_kegg_id):
        
        # Call the function

        # Run the function with test data
        contig_to_kegg_counter, contig_to_genes = manage_protein_alignement(
            faa_file=Path(faa_file),
            contigs_fasta=Path("contigs_fasta"),
            contig_to_length=contig_to_length,
            contigs_in_bins=set(),
            diamond_result_file=Path("diamond_result_file"),
            checkm2_db=None,
            threads=1,
            resume=True,
            low_mem=False
        )

    # Assertions to check the function output or file existence
    assert isinstance(contig_to_genes, dict)
    assert isinstance(contig_to_kegg_counter, dict)
    assert len(contig_to_genes) == 3


def test_manage_protein_alignement_not_resume(tmpdir, tmp_path):
    # Create temporary directories and files for testing

    faa_file = tmp_path / "proteins.faa"
    faa_file_content = (
    ">contig1_1\nACGT\n>contig2_1\nTGCA\n>contig2_2\nAAAA\n>contig3_1\nCCCC\n"
    )

    contig_to_length={"contig1":40, "contig2":80, "contig3":20}

    faa_file.write_text(faa_file_content)


    contigs_fasta = os.path.join(str(tmpdir), "contigs.fasta")
    diamond_result_file = os.path.join(str(tmpdir), "diamond_results.tsv")

    contig_to_kegg_id = {
        "contig1": Counter({"K12345": 1, "K67890": 1}),
        "contig2": Counter({"K23456": 1})
    }


    with patch("binette.diamond.get_contig_to_kegg_id", return_value=contig_to_kegg_id), \
                             patch("binette.diamond.run", return_value=None):
        
        # Call the function

        contig_to_kegg_counter, contig_to_genes = manage_protein_alignement(
            faa_file=Path(faa_file),
            contigs_fasta=Path(contigs_fasta),
            contig_to_length=contig_to_length,
            contigs_in_bins=set(),
            diamond_result_file=Path(diamond_result_file),
            checkm2_db=None,
            threads=1,
            resume=True,
            low_mem=False
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
    fasta_file_content = (
    ">contig1\nACGT\n>contig2\nTGCA\n>contig3\nAAAA\n>contig4\nCCCC\n>contig5\nCGTCGCT\n"
    )
    fasta_file.write_text(fasta_file_content)

    # Call the function and capture the return values
    bin_set_name_to_bins, original_bins, contigs_in_bins, contig_to_length = parse_input_files(None, [bin_set1, bin_set2], fasta_file)


    # # Perform assertions on the returned values
    assert isinstance(bin_set_name_to_bins, dict)
    assert isinstance(original_bins, set)
    assert isinstance(contigs_in_bins, set)
    assert isinstance(contig_to_length, dict)


    assert set(bin_set_name_to_bins) == {'bin_set1', "bin_set2"}
    assert len(original_bins) == 4
    assert contigs_in_bins == {"contig1","contig2", "contig3","contig4"}
    assert len(contig_to_length) == 4

def test_parse_input_files_with_contig2bin_tables_with_unknown_contig(tmp_path):

    bin_set3 = tmp_path / "bin_set3.tsv"
    bin_set3.write_text("contig3\tbin3A\ncontig44\ttbin3B\n")
    fasta_file = tmp_path / "assembly.fasta"
    fasta_file_content = (
    ">contig1\nACGT\n>contig2\nTGCA\n>contig3\nAAAA\n>contig4\nCCCC\n>contig5\nCGTCGCT\n"
    )
    fasta_file.write_text(fasta_file_content)

    with pytest.raises(ValueError):
        parse_input_files(None, [bin_set3], fasta_file)


def test_parse_input_files_bin_dirs(create_temp_bin_directories, tmp_path):

    bin_dirs = [Path(d) for d in create_temp_bin_directories.values()]

    contig2bin_tables = []

    # Create temporary directories and files for testing

    fasta_file = tmp_path / "assembly.fasta"
    fasta_file_content = (
    ">contig1\nACGT\n>contig2\nTGCA\n>contig3\nAAAA\n>contig4\nCCCC\n>contig5\nCGTCGCT\n"
    )
    fasta_file.write_text(fasta_file_content)

    # Call the function and capture the return values
    bin_set_name_to_bins, original_bins, contigs_in_bins, contig_to_length = parse_input_files(bin_dirs, contig2bin_tables, fasta_file)

    # # Perform assertions on the returned values
    assert isinstance(bin_set_name_to_bins, dict)
    assert isinstance(original_bins, set)
    assert isinstance(contigs_in_bins, set)
    assert isinstance(contig_to_length, dict)


    assert set(bin_set_name_to_bins) == {'set1', 'set2'}
    assert len(original_bins) == 3
    assert contigs_in_bins == {"contig1","contig2", "contig3","contig4","contig5",}
    assert len(contig_to_length) == 5


def test_argument_used_once():
    # Test UniqueStore class 
    parser = ArgumentParser(description='Test parser')
    parser.add_argument('--example', action=UniqueStore, help='Example argument')
    args = parser.parse_args(['--example', 'value'])
    assert args.example == 'value'

def test_argument_used_multiple_times():
    # Test UniqueStore class 
    parser = ArgumentParser(description='Test parser')
    parser.add_argument('--example', action=UniqueStore, help='Example argument')
    with pytest.raises(SystemExit):
        parser.parse_args(['--example', 'value', '--example', 'value2'])


def test_parse_arguments_required_arguments():
    # Test when only required arguments are provided
    args = parse_arguments(["-d", "folder1", "folder2", "-c", "contigs.fasta"])
    assert args.bin_dirs == [Path("folder1"), Path("folder2")]
    assert args.contigs == Path("contigs.fasta")

def test_parse_arguments_optional_arguments():
    # Test when required and optional arguments are provided
    args = parse_arguments(["-d", "folder1", "folder2", "-c", "contigs.fasta", "--threads", "4", "--outdir", "output"])
    assert args.bin_dirs == [Path("folder1"), Path("folder2")]
    assert args.contigs == Path("contigs.fasta")
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
    with patch('binette.contig_manager.parse_fasta_file') as mock_parse_fasta_file, \
         patch('binette.cds.predict') as mock_predict, \
         patch('binette.diamond.get_checkm2_db') as mock_get_checkm2_db, \
         patch('binette.diamond.run') as mock_diamond_run, \
         patch('binette.diamond.get_contig_to_kegg_id') as mock_diamond_get_contig_to_kegg_id:
        
        # Set the return value of the mocked functions
        mock_parse_fasta_file.return_value = [MagicMock(name="contig1")]
        mock_predict.return_value = {"contig1": ["gene1"]}
        
        # Call the function
        contig_to_kegg_counter, contig_to_genes = manage_protein_alignement(
            faa_file, contigs_fasta, contig_to_length, contigs_in_bins,
            diamond_result_file, checkm2_db, threads, resume, low_mem
        )
        
        # Assertions to check if functions were called
        mock_parse_fasta_file.assert_called_once_with(contigs_fasta.as_posix())
        mock_predict.assert_called_once()
        mock_diamond_get_contig_to_kegg_id.assert_called_once()
        mock_diamond_run.assert_called_once_with(
            faa_file.as_posix(), diamond_result_file.as_posix(), checkm2_db.as_posix(), f"{os.path.splitext(diamond_result_file.as_posix())[0]}.log", threads, low_mem=low_mem
        )

def test_main_resume_when_not_possible(monkeypatch):
    # Define or mock the necessary inputs/arguments

    # Mock sys.argv to use test_args
    test_args = [
        "-d", "bin_dir1", "bin_dir2",
        "-c", "contigs.fasta",
        # ... more arguments as required ...
        "--debug",
        "--resume"
    ]
    monkeypatch.setattr(sys, 'argv', ['your_script.py'] + test_args)

    # Call the main function
    with pytest.raises(FileNotFoundError):
        main()

def test_main(monkeypatch):
    # Define or mock the necessary inputs/arguments

    # Mock sys.argv to use test_args
    test_args = [
        "-d", "bin_dir1", "bin_dir2",
        "-c", "contigs.fasta",
        # ... more arguments as required ...
        "--debug"
    ]
    monkeypatch.setattr(sys, 'argv', ['your_script.py'] + test_args)

    # Mock the necessary functions
    with patch('binette.main.parse_input_files') as mock_parse_input_files, \
         patch('binette.main.manage_protein_alignement') as mock_manage_protein_alignement, \
         patch('binette.contig_manager.apply_contig_index') as mock_apply_contig_index, \
         patch('binette.bin_manager.rename_bin_contigs') as mock_rename_bin_contigs, \
         patch('binette.bin_manager.create_intermediate_bins') as mock_create_intermediate_bins, \
         patch('binette.bin_quality.add_bin_metrics') as mock_add_bin_metrics, \
         patch('binette.main.log_selected_bin_info') as mock_log_selected_bin_info, \
         patch('binette.contig_manager.make_contig_index') as mock_make_contig_index, \
         patch('binette.io_manager.write_original_bin_metrics') as mock_write_original_bin_metrics, \
         patch('binette.main.select_bins_and_write_them') as mock_select_bins_and_write_them:
        
        # Set return values for mocked functions if needed
        mock_parse_input_files.return_value = (None, None, None, None)
        mock_manage_protein_alignement.return_value = ({"contig1": 1}, {"contig1": ["gene1"]})
        mock_make_contig_index.return_value = ({}, {})
        mock_apply_contig_index.return_value = MagicMock()
        mock_rename_bin_contigs.return_value = MagicMock()
        mock_create_intermediate_bins.return_value = MagicMock()
        mock_add_bin_metrics.return_value = MagicMock()
        mock_log_selected_bin_info.return_value = MagicMock()
        
        
        main()


        # Add assertions to ensure the mocks were called as expected
        mock_parse_input_files.assert_called_once()
        mock_manage_protein_alignement.assert_called_once()
        mock_rename_bin_contigs.assert_called_once()
        mock_create_intermediate_bins.assert_called_once()
        
        mock_log_selected_bin_info.assert_called_once()
        mock_select_bins_and_write_them.assert_called_once()
        mock_write_original_bin_metrics.assert_called_once()

        assert mock_apply_contig_index.call_count == 3
        assert mock_add_bin_metrics.call_count == 2