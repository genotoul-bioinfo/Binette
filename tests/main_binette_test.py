
import pytest
import logging
from binette.binette import log_selected_bin_info, select_bins_and_write_them, manage_protein_alignement, parse_input_files, parse_arguments, init_logging, main
from binette.bin_manager import Bin
from binette import diamond
import os
import sys
from unittest.mock import patch

from collections import Counter
from bin_manager_test import create_temp_bin_directories, create_temp_bin_files

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
        set(bins), str(contigs_fasta), final_bin_report, min_completeness=60, index_to_contig=index_to_contig, outdir=str(outdir), debug=False
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
            faa_file=str(faa_file),
            contigs_fasta="contigs_fasta",
            contig_to_length=contig_to_length,
            contigs_in_bins={},
            diamond_result_file="diamond_result_file",
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
            faa_file=str(faa_file),
            contigs_fasta=contigs_fasta,
            contig_to_length=contig_to_length,
            contigs_in_bins={},
            diamond_result_file=diamond_result_file,
            checkm2_db=None,
            threads=1,
            resume=True,
            low_mem=False
        )

    # Assertions to check the function output or file existence
    assert isinstance(contig_to_genes, dict)
    assert isinstance(contig_to_kegg_counter, dict)
    assert len(contig_to_genes) == 3


def test_parse_input_files_bin_dirs(create_temp_bin_directories, tmp_path):

    set_name_to_bin_dir = create_temp_bin_directories
    bin_dirs = list(create_temp_bin_directories.values())

    contig2bin_tables = []

    # Create temporary directories and files for testing

    fasta_file = tmp_path / "assembly.fasta"
    fasta_file_content = (
    ">contig1\nACGT\n>contig2\nTGCA\n>contig3\nAAAA\n>contig4\nCCCC\n>contig5\nCGTCGCT\n"
    )
    fasta_file.write_text(fasta_file_content)

    # Call the function and capture the return values
    bin_set_name_to_bins, original_bins, contigs_in_bins, contig_to_length = parse_input_files(bin_dirs, contig2bin_tables, str(fasta_file))


    # # Perform assertions on the returned values
    assert isinstance(bin_set_name_to_bins, dict)
    assert isinstance(original_bins, set)
    assert isinstance(contigs_in_bins, set)
    assert isinstance(contig_to_length, dict)


    assert set(bin_set_name_to_bins) == {'1', "2"}
    assert len(original_bins) == 3
    assert contigs_in_bins == {"contig1","contig2", "contig3","contig4","contig5",}
    assert len(contig_to_length) == 5


def test_parse_input_files_bin_dirs(create_temp_bin_directories, tmp_path):

    set_name_to_bin_dir = create_temp_bin_directories
    bin_dirs = list(create_temp_bin_directories.values())

    contig2bin_tables = []

    # Create temporary directories and files for testing

    fasta_file = tmp_path / "assembly.fasta"
    fasta_file_content = (
    ">contig1\nACGT\n>contig2\nTGCA\n>contig3\nAAAA\n>contig4\nCCCC\n>contig5\nCGTCGCT\n"
    )
    fasta_file.write_text(fasta_file_content)

    # Call the function and capture the return values
    bin_set_name_to_bins, original_bins, contigs_in_bins, contig_to_length = parse_input_files(bin_dirs, contig2bin_tables, str(fasta_file))

    # # Perform assertions on the returned values
    assert isinstance(bin_set_name_to_bins, dict)
    assert isinstance(original_bins, set)
    assert isinstance(contigs_in_bins, set)
    assert isinstance(contig_to_length, dict)


    assert set(bin_set_name_to_bins) == {'1', "2"}
    assert len(original_bins) == 3
    assert contigs_in_bins == {"contig1","contig2", "contig3","contig4","contig5",}
    assert len(contig_to_length) == 5



def test_parse_arguments_required_arguments():
    # Test when only required arguments are provided
    args = parse_arguments(["-d", "folder1", "folder2", "-c", "contigs.fasta"])
    assert args.bin_dirs == ["folder1", "folder2"]
    assert args.contigs == "contigs.fasta"

def test_parse_arguments_optional_arguments():
    # Test when required and optional arguments are provided
    args = parse_arguments(["-d", "folder1", "folder2", "-c", "contigs.fasta", "--threads", "4", "--outdir", "output"])
    assert args.bin_dirs == ["folder1", "folder2"]
    assert args.contigs == "contigs.fasta"
    assert args.threads == 4
    assert args.outdir == "output"

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



def test_main_functionality_resume_when_not_possible(monkeypatch):
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

    # You may also need to mock certain functions to avoid actual file operations or to simulate their behavior
    # For example, mock the functions parse_input_files, manage_protein_alignement, select_bins_and_write_them, etc.

    # Call the main function
    with pytest.raises(FileNotFoundError) as e_info:
        main()


# def test_main_functionality_(monkeypatch, create_temp_bin_directories, tmp_path):
#     # Define or mock the necessary inputs/arguments

#     bin_dirs = list(create_temp_bin_directories.values())
#     fasta_file = tmp_path / "assembly.fasta"
#     fasta_file_content = (
#     ">contig1\nACGT\n>contig2\nTGCA\n>contig3\nAAAA\n>contig4\nCCCC\n>contig5\nCGTCGCT\n"
#     )
#     fasta_file.write_text(fasta_file_content)

#     # Mock sys.argv to use test_args
#     test_args = [
#         "-d"] + bin_dirs + ["-c", str(fasta_file),
#         # ... more arguments as required ...
#         "--debug"
#     ]

#     monkeypatch.setattr(sys, 'argv', ['your_script.py'] + test_args)

#     # You may also need to mock certain functions to avoid actual file operations or to simulate their behavior
#     # For example, mock the functions parse_input_files, manage_protein_alignement, select_bins_and_write_them, etc.

#     # Call the main function
#     with patch("binette.diamond.get_contig_to_kegg_id", return_value=contig_to_kegg_id), \
#                             patch("binette.diamond.run", return_value=None):
#         result = main()
#         assert result == 0
        