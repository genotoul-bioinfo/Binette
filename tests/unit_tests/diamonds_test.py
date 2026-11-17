import shutil
import subprocess
from collections import Counter
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from binette import diamond


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
    with pytest.raises(diamond.DiamondError) as exc_info:
        diamond.get_checkm2_db()

    assert exc_info.value.exit_code == 1


def test_get_checkm2_db_with_success(monkeypatch):
    def mock_subprocess_run(*args, **kwargs):
        # Simulating the behavior of checkm2 command
        if (
            args[0][0] == "checkm2"
            and args[0][1] == "database"
            and args[0][2] == "--current"
        ):
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
        if (
            args[0][0] == "checkm2"
            and args[0][1] == "database"
            and args[0][2] == "--current"
        ):
            return CompletedProcess(2, "")

    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)
    monkeypatch.setattr(shutil, "which", mock_shutil_which)

    # Call the function
    with pytest.raises(diamond.DiamondError) as exc_info:
        diamond.get_checkm2_db()

    assert exc_info.value.exit_code == 1


def test_get_checkm2_db_wrong_path_format(monkeypatch):
    def mock_subprocess_run(*args, **kwargs):
        # Simulating the behavior of checkm2 command
        if (
            args[0][0] == "checkm2"
            and args[0][1] == "database"
            and args[0][2] == "--current"
        ):
            return CompletedProcess(
                0,
                "UNEXPECTED PATH FORMAT RETURNED BY CHECKM2: /mocked/path/to/checkm2.dmnd",
            )

    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)
    monkeypatch.setattr(shutil, "which", mock_shutil_which)

    # Call the function
    with pytest.raises(diamond.DiamondError) as exc_info:
        diamond.get_checkm2_db()

    assert exc_info.value.exit_code == 1


def test_check_tool_exists_tool_found(monkeypatch):
    # Mocking shutil.which
    def mock_shutil_which(*args, **kwargs):
        return "path/to/tool"

    monkeypatch.setattr(shutil, "which", mock_shutil_which)

    # Call the function
    diamond.check_tool_exists("existing_tool")


def test_check_tool_exists_tool_not_found(monkeypatch):
    # Mocking shutil.which
    def mock_shutil_which(*args, **kwargs):
        return None

    monkeypatch.setattr(shutil, "which", mock_shutil_which)

    with pytest.raises(diamond.DiamondError) as exc_info:
        diamond.check_tool_exists("non_existing_tool")

    assert exc_info.value.exit_code == 1
    assert "non_existing_tool" in str(exc_info.value)


def test_run_diamond_tool_found(monkeypatch, tmp_path):
    log_path = tmp_path / "log.txt"
    output_file = tmp_path / "output.txt"
    partial_output_file = tmp_path / "partial_output.txt"

    # Mocking subprocess.run
    def mock_subprocess_run(*args, **kwargs):
        class CompletedProcess:
            def __init__(self, returncode):
                self.returncode = returncode

        # DIAMOND writes to a temporary file that is moved on success
        expected_cmd = (
            "diamond blastp --outfmt 6 --max-target-seqs 1 --query test.faa "
            f"-o {partial_output_file.as_posix()} --threads 1 --db db --compress 1 "
            "--query-cover 80 --subject-cover 80 --id 30 --evalue 1e-05 "
            f"--block-size 2 2> {log_path.as_posix()}"
        )
        assert args[0] == expected_cmd
        partial_output_file.write_text("diamond results")
        return CompletedProcess(0)

    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)

    # Call the function
    diamond.run(
        "test.faa",
        output_file.as_posix(),
        "db",
        log_path.as_posix(),
        threads=1,
        query_cover=80,
        subject_cover=80,
        percent_id=30,
        evalue=1e-05,
        low_mem=False,
    )

    assert output_file.read_text() == "diamond results"
    assert not partial_output_file.exists()


def test_run_diamond_fails(monkeypatch, tmp_path):
    log_file = tmp_path / "log.txt"
    output_file = tmp_path / "output.txt"

    def mock_subprocess_run(*args, **kwargs):
        class CompletedProcess:
            returncode = 1

        # Simulate a partial output left by a failing DIAMOND run
        (tmp_path / "partial_output.txt").write_text("truncated")
        return CompletedProcess()

    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)
    monkeypatch.setattr(diamond, "check_tool_exists", lambda tool: None)

    with pytest.raises(diamond.DiamondError) as exc_info:
        diamond.run(
            "test.faa",
            output_file.as_posix(),
            "db",
            log_file.as_posix(),
            threads=1,
        )

    assert exc_info.value.exit_code == 1
    # A failed run never produces the final result file
    assert not output_file.exists()


def test_get_contig_to_kegg_id():
    # Mock input data
    diamond_result_file = "dummy_diamond_results.txt"

    # Mocked dataframe representing the data read from the Diamond result file
    mocked_data = {
        "ProteinID": [
            "contig1_protein1",
            "contig1_protein2",
            "contig2_protein1",
            "contig2_protein2",
        ],
        "annotation": [
            "protein1_annotation~K12345",
            "protein2_annotation~K67890",
            "protein3_annotation~K23456",
            "protein4_annotation~K66666",
        ],
    }
    mocked_df = pd.DataFrame(mocked_data)

    # Mocked return values for keggData.KeggCalculator() and KeggCalc.return_default_values_from_category()
    class MockedKeggCalculator:
        def return_default_values_from_category(self, category):
            return {"K12345": 2, "K67890": 1, "K23456": 3}

    mocked_kegg_calculator = MockedKeggCalculator()

    # Mocking relevant functions and classes used within the function
    with (
        patch("pandas.read_csv", return_value=mocked_df),
        patch("checkm2.keggData.KeggCalculator", return_value=mocked_kegg_calculator),
    ):
        # Call the function
        result = diamond.get_contig_to_kegg_id(diamond_result_file)

    # K66666 is not in return_default_values_from_category so it won't be kept in result kegg

    # Define the expected output based on the mocked data
    expected_result = {
        "contig1": Counter({"K12345": 1, "K67890": 1}),
        "contig2": Counter({"K23456": 1}),
    }

    # Check if the function output matches the expected result
    assert result == expected_result


def test_get_contig_to_kegg_id_empty_file():
    """Test that get_contig_to_kegg_id raises an error with exit code 3 when DIAMOND result file is empty."""
    diamond_result_file = "empty_diamond_results.txt"

    # Mock empty dataframe
    empty_df = pd.DataFrame()

    with patch("pandas.read_csv", return_value=empty_df):
        with pytest.raises(diamond.DiamondError) as exc_info:
            diamond.get_contig_to_kegg_id(diamond_result_file)

    assert exc_info.value.exit_code == 3


def test_resolve_checkm2_db(tmp_path):
    db = tmp_path / "db.dmnd"
    db.write_text("db")
    assert diamond.resolve_checkm2_db(db) == db

    with pytest.raises(FileNotFoundError):
        diamond.resolve_checkm2_db(tmp_path / "missing.dmnd")

    with patch("binette.diamond.get_checkm2_db", return_value="/some/db.dmnd"):
        assert diamond.resolve_checkm2_db(None) == Path("/some/db.dmnd")


def test_result_parameters_include_block_size():
    default = diamond.result_parameters(low_mem=False)
    low_mem = diamond.result_parameters(low_mem=True)

    assert default["block_size"] == 2
    assert low_mem["block_size"] == 0.5
    assert {k: v for k, v in default.items() if k != "block_size"} == {
        k: v for k, v in low_mem.items() if k != "block_size"
    }
    assert {"query_cover", "subject_cover", "percent_id", "evalue"} <= set(default)
