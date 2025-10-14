import logging
from pathlib import Path

import pytest
from typer.testing import CliRunner

from binette.main import app

runner = CliRunner()
logger = logging.getLogger(__name__)


def test_help_app():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0


def test_version_app():
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert "Binette" in result.output


def test_no_bin_input(tmp_path):
    input = tmp_path / "file.txt"
    input.write_text("This is not a fasta file")
    result = runner.invoke(app, ["-c", input])

    assert result.exit_code == 2


def test_wrong_input(tmp_path):
    input = tmp_path / "file.txt"
    input.write_text("This is not a fasta file")
    result = runner.invoke(app, ["-c", input, "--bin_dirs", tmp_path])

    print(result.output)
    print(result.stderr)
    print(result.exception)
    assert f"{input.as_posix()} is not fasta or fastq sequence file" in str(
        result.exception
    )
    assert result.exit_code == 1


def test_wrong_input_both_input_types(tmp_path):
    input = tmp_path / "file.txt"
    input.write_text("This is not a fasta file")
    result = runner.invoke(
        app,
        [
            "-c",
            input,
            "--bin_dirs",
            tmp_path,
            "--contig2bin_tables",
            tmp_path,
        ],
    )

    print("STDOUT:", result.output)
    print("STDERR:", result.stderr)

    print("EXCEPTION:", result.exception)

    assert result.exit_code == 2


def test_quiet_and_verbose_flag(tmp_path):
    input = tmp_path / "contig.fna"
    input.write_text(">contig1\nATGC\n>contig2\nATGC")
    result = runner.invoke(app, ["-c", input, "--bin_dirs", tmp_path, "-q", "-v"])

    print("STDOUT:", result.output)
    print("STDERR:", result.stderr)

    print("EXCEPTION:", result.exception)
    assert result.exit_code == 2


def compare_results(result_file: Path, expected_file: Path):
    content_actual = result_file.read_text()
    content_expected = expected_file.read_text()
    print(content_actual)
    assert content_actual == content_expected, (
        f"Content mismatch for {result_file}. Expected content from {expected_file}."
    )


@pytest.mark.requires_test_data
def test_bin_dir_input(test_data_path: Path, tmp_path):
    binning_results_dir = test_data_path / "binning_results"

    cmd_args = [
        "-d",
        binning_results_dir / "A/",
        "-d",
        binning_results_dir / "B/",
        "-d",
        binning_results_dir / "C/",
        "--contigs",
        test_data_path / "all_contigs.fna",
        "--checkm2_db",
        test_data_path / "checkm2_tiny_db/checkm2_tiny_db.dmnd",
        "-o",
        tmp_path / "test_results_from_dirs",
    ]
    cmd_args = [arg.as_posix() if isinstance(arg, Path) else arg for arg in cmd_args]

    result = runner.invoke(app, cmd_args)
    print(result.output)
    print(result.stderr)
    assert result.exit_code == 0

    result_table = (
        tmp_path / "test_results_from_dirs" / "final_bins_quality_reports.tsv"
    )
    expected_table = Path("tests/expected_results/final_bins_quality_reports.tsv")
    compare_results(result_table, expected_table)


@pytest.mark.requires_test_data
def test_bin_tables_input_and_resume(test_data_path: Path, tmp_path):
    binning_results_dir = test_data_path / "binning_results"

    cmd_args = [
        "-b",
        binning_results_dir / "A.binning",
        "-b",
        binning_results_dir / "B.binning",
        "-b",
        binning_results_dir / "C.binning",
        "--contigs",
        test_data_path / "all_contigs.fna",
        "--checkm2_db",
        test_data_path / "checkm2_tiny_db/checkm2_tiny_db.dmnd",
        "-o",
        tmp_path / "test_results_from_dirs",
        "-q",
    ]
    cmd_args = [arg.as_posix() if isinstance(arg, Path) else arg for arg in cmd_args]

    result = runner.invoke(app, cmd_args)
    print(result.output)
    print(result.stderr)

    assert result.exit_code == 0

    result_table = (
        tmp_path / "test_results_from_dirs" / "final_bins_quality_reports.tsv"
    )
    expected_table = Path("tests/expected_results/final_bins_quality_reports.tsv")
    compare_results(result_table, expected_table)

    cmd_args.append("--resume")
    result = runner.invoke(app, cmd_args)
    print(result.output)
    print(result.stderr)
    assert result.exit_code == 0

    expected_table = Path(
        "tests/expected_results/final_bins_quality_reports_from_proteins_input.tsv"
    )
    compare_results(result_table, expected_table)


@pytest.mark.requires_test_data
def test_bin_tables_input_and_protein_input(test_data_path: Path, tmp_path):
    binning_results_dir = test_data_path / "binning_results"

    cmd_args = [
        "-b",
        binning_results_dir / "A.binning",
        "-b",
        binning_results_dir / "B.binning",
        "-b",
        binning_results_dir / "C.binning",
        "--contigs",
        test_data_path / "all_contigs.fna",
        "--checkm2_db",
        test_data_path / "checkm2_tiny_db/checkm2_tiny_db.dmnd",
        "--proteins",
        test_data_path / "proteins.faa",
        "-o",
        tmp_path / "test_results_from_dirs",
        "-v",
        "--debug",
    ]
    cmd_args = [arg.as_posix() if isinstance(arg, Path) else arg for arg in cmd_args]

    result = runner.invoke(app, cmd_args)
    print(result.output)
    print(result.stderr)

    assert result.exit_code == 0

    result_table = (
        tmp_path / "test_results_from_dirs" / "final_bins_quality_reports.tsv"
    )
    expected_table = Path(
        "tests/expected_results/final_bins_quality_reports_from_proteins_input.tsv"
    )
    compare_results(result_table, expected_table)
