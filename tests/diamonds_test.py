import re
import subprocess
import shutil
import sys
import logging

from unittest.mock import patch, MagicMock

import subprocess
import sys
import logging
import pytest

from binette import diamond

import pandas as pd
from collections import Counter

class CompletedProcess:
    def __init__(self, returncode, stderr):
        self.returncode = returncode
        self.stderr = stderr


def mock_shutil_which(*args, **kwargs):
    if args[0] == "checkm2":
        return "checkm2"

def test_get_checkm2_db_no_checkm2(monkeypatch):
    # Mocking shutil.which
    def mock_shutil_which_none(*args, **kwargs):
        return None

    monkeypatch.setattr(shutil, "which", mock_shutil_which_none)

    # Call the function
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        diamond.get_checkm2_db()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1

def test_get_checkm2_db_with_success(monkeypatch):

    def mock_subprocess_run(*args, **kwargs):

        # Simulating the behavior of checkm2 command
        if args[0][0] == "checkm2" and args[0][1] == "database" and args[0][2] == "--current":
            return CompletedProcess(0, "INFO: /mocked/path/to/checkm2.dmnd")
            
    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)
    monkeypatch.setattr(shutil, "which", mock_shutil_which)

    # Call the function
    result = diamond.get_checkm2_db()

    expected_path = "/mocked/path/to/checkm2.dmnd"
    assert result == expected_path

def test_get_checkm2_db_checkm2_exit_error(monkeypatch):

    def mock_subprocess_run(*args, **kwargs):

        # Simulating the behavior of checkm2 command
        if args[0][0] == "checkm2" and args[0][1] == "database" and args[0][2] == "--current":
            return CompletedProcess(2, "")
            

    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)
    monkeypatch.setattr(shutil, "which", mock_shutil_which)

    # Call the function
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        diamond.get_checkm2_db()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1

def test_get_checkm2_db_wrong_path_format(monkeypatch):

    def mock_subprocess_run(*args, **kwargs):

        # Simulating the behavior of checkm2 command
        if args[0][0] == "checkm2" and args[0][1] == "database" and args[0][2] == "--current":
            return CompletedProcess(0, "UNEXPECTED PATH FORMAT RETURNED BY CHECKM2: /mocked/path/to/checkm2.dmnd")
            
    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)
    monkeypatch.setattr(shutil, "which", mock_shutil_which)

    # Call the function
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        diamond.get_checkm2_db()

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1



def test_check_tool_exists_tool_found(monkeypatch):
    # Mocking shutil.which
    def mock_shutil_which(*args, **kwargs):
        return "path/to/tool"

    monkeypatch.setattr(shutil, "which", mock_shutil_which)

    # Call the function
    try:
        diamond.check_tool_exists("existing_tool")
    except FileNotFoundError:
        pytest.fail("check_tool_exists raised FileNotFoundError unexpectedly.")


def test_check_tool_exists_tool_not_found(monkeypatch):
    # Mocking shutil.which
    def mock_shutil_which(*args, **kwargs):
        return None

    monkeypatch.setattr(shutil, "which", mock_shutil_which)

    # Call the function and expect FileNotFoundError
    with pytest.raises(FileNotFoundError):
        diamond.check_tool_exists("non_existing_tool")
        




def test_run_diamond_tool_found(monkeypatch):

    monkeypatch.setattr(sys, "exit", lambda x: None)  # Patch sys.exit to avoid test interruption

    # Mocking subprocess.run
    def mock_subprocess_run(*args, **kwargs):
        class CompletedProcess:
            def __init__(self, returncode):
                self.returncode = returncode

        # Simulating successful run of diamond command
        if args[0] == "diamond blastp --outfmt 6 --max-target-seqs 1 --query test.faa -o output.txt --threads 1 --db db --query-cover 80 --subject-cover 80 --id 30 --evalue 1e-05 --block-size 2 2> log.txt":
            return CompletedProcess(0)

    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)
    monkeypatch.setattr(logging, "error", lambda x: None)  # Avoid logging during test

    # Call the function
    diamond.run(
        "test.faa", "output.txt", "db", "log.txt", threads=1, query_cover=80, subject_cover=80, percent_id=30,
        evalue=1e-05, low_mem=False
    )


def test_run_diamond_tool_not_found(monkeypatch):
    # Mocking check_tool_exists to simulate tool not found scenario
    def mock_check_tool_exists(*args, **kwargs):
        raise FileNotFoundError

    monkeypatch.setattr(logging, "error", lambda x: None)  # Avoid logging during test

    # Call the function and expect it to raise FileNotFoundError
    with patch('sys.exit') as mock_exit:
        diamond.run(
            "test.faa", "output.txt", "db", "log.txt", threads=1, query_cover=80, subject_cover=80, percent_id=30,
            evalue=1e-05, low_mem=False
        )

    mock_exit.assert_called_once_with(1)


def test_get_contig_to_kegg_id():
    # Mock input data
    diamond_result_file = "dummy_diamond_results.txt"

    # Mocked dataframe representing the data read from the Diamond result file
    mocked_data = {
        "ProteinID": ["contig1_protein1", "contig1_protein2", "contig2_protein1", "contig2_protein2"],
        "annotation": ["protein1_annotation~K12345", "protein2_annotation~K67890", "protein3_annotation~K23456", "protein4_annotation~K66666"]
    }
    mocked_df = pd.DataFrame(mocked_data)

    # Mocked return values for keggData.KeggCalculator() and KeggCalc.return_default_values_from_category()
    class MockedKeggCalculator:
        def return_default_values_from_category(self, category):
            return {"K12345": 2, "K67890": 1, "K23456":3}

    mocked_kegg_calculator =  MockedKeggCalculator()

    # Mocking relevant functions and classes used within the function
    with patch("pandas.read_csv", return_value=mocked_df), \
         patch("checkm2.keggData.KeggCalculator", return_value=mocked_kegg_calculator):
        
        # Call the function
        result = diamond.get_contig_to_kegg_id(diamond_result_file)

    # K66666 is not in return_default_values_from_category so it won't be kept in result kegg

    # Define the expected output based on the mocked data
    expected_result = {
        "contig1": Counter({"K12345": 1, "K67890": 1}),
        "contig2": Counter({"K23456": 1})
    }

    # Check if the function output matches the expected result
    assert result == expected_result


# Additional tests can be added to cover more edge cases and scenarios.
