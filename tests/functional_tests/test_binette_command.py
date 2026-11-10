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


def bin_tables_cmd(test_data_path: Path, outdir: Path, bin_sets=("A", "B", "C")):
    binning_results_dir = test_data_path / "binning_results"
    cmd_args = []
    for bin_set in bin_sets:
        cmd_args += ["-b", binning_results_dir / f"{bin_set}.binning"]
    cmd_args += [
        "--contigs",
        test_data_path / "all_contigs.fna",
        "--checkm2_db",
        test_data_path / "checkm2_tiny_db/checkm2_tiny_db.dmnd",
        "-o",
        outdir,
    ]
    return [arg.as_posix() if isinstance(arg, Path) else arg for arg in cmd_args]


def run_binette(cmd_args, caplog):
    caplog.clear()
    result = runner.invoke(app, cmd_args)
    print(result.output)
    print(result.stderr)
    print(caplog.text)
    return result


@pytest.mark.requires_test_data
def test_bin_tables_input_and_resume(test_data_path: Path, tmp_path, caplog):
    caplog.set_level(logging.INFO, logger="binette")
    outdir = tmp_path / "test_results"
    cmd_args = bin_tables_cmd(test_data_path, outdir)

    result_table = outdir / "final_bins_quality_reports.tsv"
    expected_table = Path("tests/expected_results/final_bins_quality_reports.tsv")
    diamond_file = outdir / "temporary_files" / "diamond_result.tsv.gz"
    faa_file = outdir / "temporary_files" / "assembly_proteins.faa.gz"

    # Fresh run: resume is on by default and nothing can be resumed
    result = run_binette(cmd_args, caplog)
    assert result.exit_code == 0
    assert "No intermediate results from a previous run" in caplog.text
    assert "Running diamond" in caplog.text
    compare_results(result_table, expected_table)
    faa_mtime = faa_file.stat().st_mtime_ns
    diamond_mtime = diamond_file.stat().st_mtime_ns
    report_mtime = result_table.stat().st_mtime_ns

    # Same inputs and parameters: nothing to do
    result = run_binette(cmd_args, caplog)
    assert result.exit_code == 0
    assert "Binette already completed" in caplog.text
    assert result_table.stat().st_mtime_ns == report_mtime

    # Explicit --resume behaves the same
    result = run_binette(cmd_args + ["--resume"], caplog)
    assert result.exit_code == 0
    assert "Binette already completed" in caplog.text

    # --no-write-fasta-bins does not change results: nothing to do
    final_bins_dir = outdir / "final_bins"
    bins_fasta = {p.name: p.read_text() for p in final_bins_dir.glob("*.fa")}
    assert bins_fasta
    result = run_binette(cmd_args + ["--no-write-fasta-bins"], caplog)
    assert result.exit_code == 0
    assert "Binette already completed" in caplog.text
    assert "Creating intermediate bins" not in caplog.text
    assert {p.name: p.read_text() for p in final_bins_dir.glob("*.fa")} == bins_fasta

    # Missing bin FASTA file: rewritten without recomputing bin selection
    removed = sorted(final_bins_dir.glob("*.fa"))[0]
    removed.unlink()
    result = run_binette(cmd_args, caplog)
    assert result.exit_code == 0
    assert "Binette already completed" in caplog.text
    assert "Creating intermediate bins" not in caplog.text
    assert {p.name: p.read_text() for p in final_bins_dir.glob("*.fa")} == bins_fasta

    # Different selection parameter: proteins and DIAMOND results are reused
    result = run_binette(cmd_args + ["--prefix", "other"], caplog)
    assert result.exit_code == 0
    assert "Resuming: reusing protein sequences" in caplog.text
    assert "Resuming: reusing DIAMOND results" in caplog.text
    assert "Running diamond" not in caplog.text
    assert faa_file.stat().st_mtime_ns == faa_mtime
    assert diamond_file.stat().st_mtime_ns == diamond_mtime
    # Bins of the previous run (other prefix) are removed
    final_bins = list((outdir / "final_bins").glob("*.fa"))
    assert final_bins
    assert all(path.name.startswith("other_") for path in final_bins)

    # Back to the original parameters: same results as the fresh run
    result = run_binette(cmd_args, caplog)
    assert result.exit_code == 0
    assert "Running diamond" not in caplog.text
    compare_results(result_table, expected_table)

    # Truncated DIAMOND result and no final results (e.g. job killed while
    # DIAMOND was writing): proteins are reused and DIAMOND is run again
    diamond_file.write_bytes(diamond_file.read_bytes()[:100])
    result_table.unlink()
    result = run_binette(cmd_args, caplog)
    assert result.exit_code == 0
    assert "Cannot reuse DIAMOND results" in caplog.text
    assert "Resuming: reusing protein sequences" in caplog.text
    assert "Running diamond" in caplog.text
    compare_results(result_table, expected_table)

    # --no-resume recomputes everything
    result = run_binette(cmd_args + ["--no-resume"], caplog)
    assert result.exit_code == 0
    assert "Predicting CDS sequences" in caplog.text
    assert "Running diamond" in caplog.text
    compare_results(result_table, expected_table)


@pytest.mark.requires_test_data
def test_outdir_with_different_inputs(test_data_path: Path, tmp_path, caplog):
    caplog.set_level(logging.INFO, logger="binette")
    outdir = tmp_path / "test_results"

    result = run_binette(bin_tables_cmd(test_data_path, outdir), caplog)
    assert result.exit_code == 0

    # Different bin sets in the same output directory: refuse to run
    other_cmd = bin_tables_cmd(test_data_path, outdir, bin_sets=("A", "B"))
    result = run_binette(other_cmd, caplog)
    assert result.exit_code == 1
    assert "different inputs" in caplog.text
    assert "input bin sets" in caplog.text
    assert "--force" in caplog.text

    # --force allows it and recomputes from scratch
    result = run_binette(other_cmd + ["--force"], caplog)
    assert result.exit_code == 0
    assert "--force is set" in caplog.text
    assert "Running diamond" in caplog.text

    # First run without FASTA bins, then ask for them: only the FASTA files are written
    result = run_binette(other_cmd + ["--no-write-fasta-bins", "--no-resume"], caplog)
    assert result.exit_code == 0
    assert not list((outdir / "final_bins").glob("*.fa"))

    # The output directory now belongs to the new run
    result = run_binette(other_cmd, caplog)
    assert result.exit_code == 0
    assert "Binette already completed" in caplog.text
    assert "Creating intermediate bins" not in caplog.text
    assert list((outdir / "final_bins").glob("*.fa"))


@pytest.mark.requires_test_data
def test_outdir_from_older_binette_version(test_data_path: Path, tmp_path, caplog):
    caplog.set_level(logging.INFO, logger="binette")
    outdir = tmp_path / "test_results"
    (outdir / "temporary_files").mkdir(parents=True)
    (outdir / "temporary_files" / "diamond_result.tsv.gz").write_text("old")
    cmd_args = bin_tables_cmd(test_data_path, outdir)

    result = run_binette(cmd_args, caplog)
    assert result.exit_code == 1
    assert "cannot be verified" in caplog.text

    result = run_binette(cmd_args + ["--force"], caplog)
    assert result.exit_code == 0
    compare_results(
        outdir / "final_bins_quality_reports.tsv",
        Path("tests/expected_results/final_bins_quality_reports.tsv"),
    )


@pytest.mark.requires_test_data
def test_bin_tables_input_and_protein_input(test_data_path: Path, tmp_path, caplog):
    caplog.set_level(logging.INFO, logger="binette")
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

    # Same run again: nothing to do
    result = run_binette(cmd_args, caplog)
    assert result.exit_code == 0
    assert "Binette already completed" in caplog.text

    # Different selection parameter: the filtered proteins and DIAMOND results are reused
    result = run_binette(cmd_args + ["--min_completeness", "50"], caplog)
    assert result.exit_code == 0
    assert "Resuming: reusing protein sequences" in caplog.text
    assert "Resuming: reusing DIAMOND results" in caplog.text
