from pathlib import Path

import pytest
from pyroaring import BitMap

from binette import io_manager
from binette.bin_manager import Bin


@pytest.fixture
def bin1():
    b = Bin(contigs=BitMap({1, 3}), origin="test1", name="bin_1")
    b.score = 80
    b.N50 = 500
    return b


@pytest.fixture
def bin2():
    b = Bin(contigs=BitMap({2, 4}), origin="test2", name="bin_2")
    b.score = 75
    b.N50 = 600
    return b


def test_infer_bin_name_from_bin_inputs():
    # Mock input data
    input_bins = ["/path/to/bin1", "/path/to/bin2", "/path/to/bin3"]

    # Call the function
    result = io_manager.infer_bin_set_names_from_input_paths(
        list(map(Path, input_bins))
    )

    # Define the expected output
    expected_result = {
        "bin1": Path("/path/to/bin1"),
        "bin2": Path("/path/to/bin2"),
        "bin3": Path("/path/to/bin3"),
    }

    # Check if the output matches the expected dictionary
    assert result == expected_result


def test_infer_bin_name_from_single_path():
    # Mock input data
    input_bins = [
        "/path/to/bin1",
    ]

    # Call the function
    result = io_manager.infer_bin_set_names_from_input_paths(
        list(map(Path, input_bins))
    )

    # Define the expected output
    expected_result = {
        "/path/to/bin1": Path("/path/to/bin1"),
    }

    # Check if the output matches the expected dictionary
    assert result == expected_result


def test_infer_bin_name_from_bin_table_inputs():
    # Mock input data
    input_bins = ["/path/to/bin1.tsv", "/path/to/bin2.tsv", "/path/to/bin3.tsv"]

    # Call the function
    result = io_manager.infer_bin_set_names_from_input_paths(
        list(map(Path, input_bins))
    )

    # Define the expected output
    expected_result = {
        "bin1": Path("/path/to/bin1.tsv"),
        "bin2": Path("/path/to/bin2.tsv"),
        "bin3": Path("/path/to/bin3.tsv"),
    }

    # Check if the output matches the expected dictionary
    assert result == expected_result


def test_infer_bin_name_from_bin_table_with_different_ext():
    # Mock input data
    input_bins = ["/path/to/bin1.tsv", "/path/to/bin2.tsv", "/path/to/bin3.txt"]

    # Call the function
    result = io_manager.infer_bin_set_names_from_input_paths(
        list(map(Path, input_bins))
    )

    # Define the expected output
    expected_result = {
        "bin1.tsv": Path("/path/to/bin1.tsv"),
        "bin2.tsv": Path("/path/to/bin2.tsv"),
        "bin3.txt": Path("/path/to/bin3.txt"),
    }

    # Check if the output matches the expected dictionary
    assert result == expected_result


def test_infer_bin_name_from_bin_table_with_different_dir():
    # Mock input data
    input_bins = [
        "/path/to/bins",
        "/path2/result_bins",
        "/path2/result/bins",
    ]

    # Call the function
    result = io_manager.infer_bin_set_names_from_input_paths(
        list(map(Path, input_bins))
    )

    # Define the expected output
    expected_result = {
        "path/to/bins": Path("/path/to/bins"),
        "path2/result_bins": Path("/path2/result_bins"),
        "path2/result/bins": Path("/path2/result/bins"),
    }

    # Check if the output matches the expected dictionary
    assert result == expected_result


def test_get_paths_common_prefix_suffix():
    # Test case 1: No paths provided
    assert io_manager.get_paths_common_prefix_suffix([]) == ([], [], [])

    # # Test case 2: Single path
    assert io_manager.get_paths_common_prefix_suffix([Path("/home/user/project")]) == (
        ["/", "home", "user", "project"],
        ["/", "home", "user", "project"],
        [],
    )

    # Test case 3: Multiple paths with common prefix and suffix
    paths = [
        Path("/home/user/project/src"),
        Path("/home/user/project/docs"),
        Path("/home/user/project/tests"),
    ]
    assert io_manager.get_paths_common_prefix_suffix(paths) == (
        ["/", "home", "user", "project"],
        [],
        [],
    )

    # Test case 4: Multiple paths with no common prefix or suffix
    paths = [
        Path("/var/log/syslog"),
        Path("/usr/local/bin/python"),
        Path("/etc/nginx/nginx.conf"),
    ]
    assert io_manager.get_paths_common_prefix_suffix(paths) == (["/"], [], [])

    # Test case 5: Multiple paths with common suffix
    paths = [Path("/home/user/docs/report.txt"), Path("/home/admin/docs/report.txt")]
    assert io_manager.get_paths_common_prefix_suffix(paths) == (
        ["/", "home"],
        ["docs", "report.txt"],
        [".txt"],
    )

    # Test case 6: Paths with a deeper common prefix and suffix
    paths = [
        Path("/data/project_a/results/output.txt"),
        Path("/data/project_b/results/output.txt"),
    ]
    assert io_manager.get_paths_common_prefix_suffix(paths) == (
        ["/", "data"],
        ["results", "output.txt"],
        [".txt"],
    )

    # Test case 7: Paths with only the root as common prefix and different suffix
    paths = [Path("/project_a/output.txt"), Path("/project_b/output.txt")]
    assert io_manager.get_paths_common_prefix_suffix(paths) == (
        ["/"],
        ["output.txt"],
        [".txt"],
    )

    # Test case 8: Paths with only the root as common prefix and different suffix
    paths = [Path("/project_a/output.txt"), Path("/project_a/output.tsv")]
    assert io_manager.get_paths_common_prefix_suffix(paths) == (
        ["/", "project_a"],
        [],
        [],
    )


def test_write_bin_info(tmp_path, bin1, bin2):
    # Mock input data
    bins = [bin1, bin2]

    output_file = tmp_path / "output.tsv"

    # Call the function
    io_manager.write_bin_info(bins, output_file)

    # Check if the file was created and its content matches the expected output
    assert Path(output_file).exists()


def test_write_bin_info_add_contig(tmp_path, bin1, bin2):
    # Mock input data
    bins = [bin1, bin2]

    output_file = tmp_path / "output.tsv"

    # Call the function
    io_manager.write_bin_info(bins, output_file, add_contigs=True)

    # Check if the file was created and its content matches the expected output
    assert Path(output_file).exists()


def test_write_bins_fasta(tmp_path, bin1, bin2):
    # Mock input data
    contigs_fasta = tmp_path / "contigs.fasta"

    contigs_fasta_content = (
        ">contig1\nACGT\n>contig2\nTGCA\n>contig3\nAAAA\n>contig4\nCCCC\n"
    )
    contigs_fasta.write_text(contigs_fasta_content)

    selected_bins = [bin1, bin2]

    outdir = tmp_path / "output_bins"
    outdir.mkdir()

    # Call the function
    io_manager.write_bins_fasta(
        selected_bins,
        contigs_fasta,
        Path(outdir),
        contigs_names=["contig0", "contig1", "contig2", "contig3", "contig4"],
    )

    # Check if the files were created and their content matches the expected output
    assert (outdir / "bin_1.fa").exists()
    assert (outdir / "bin_2.fa").exists()

    with open(outdir / "bin_1.fa") as bin1_file:
        assert bin1_file.read() == ">contig1\nACGT\n>contig3\nAAAA\n"

    with open(outdir / "bin_2.fa") as bin2_file:
        assert bin2_file.read() == ">contig2\nTGCA\n>contig4\nCCCC\n"


def test_check_contig_consistency_error():
    # Mock input data
    contigs_from_assembly = ["contig1", "contig2", "contig3"]
    contigs_from_bins = ["contig2", "contig3", "contig4"]
    assembly_file = "assembly.fasta"
    elsewhere_file = "external.fasta"

    with pytest.raises(AssertionError):
        # Call the function
        io_manager.check_contig_consistency(
            contigs_from_assembly, contigs_from_bins, assembly_file, elsewhere_file
        )


def test_check_contig_consistency_no_error():
    # Mock input data
    contigs_from_assembly = ["contig1", "contig2", "contig3", "contig4"]
    contigs_from_bins = ["contig1", "contig2", "contig3"]
    assembly_file = "assembly.fasta"
    elsewhere_file = "external.fasta"

    io_manager.check_contig_consistency(
        contigs_from_assembly, contigs_from_bins, assembly_file, elsewhere_file
    )


@pytest.mark.parametrize("filename", ["coding_length.tsv", "coding_length.tsv.gz"])
def test_write_and_read_contig_coding_length(tmp_path, filename):
    output = tmp_path / filename
    contig_to_coding_length = {"contig1": 120, "contig_2": 0, "contig3": None}

    io_manager.write_contig_coding_length(contig_to_coding_length, output)

    assert output.exists()
    assert not (tmp_path / f"partial_{filename}").exists()
    assert io_manager.read_contig_coding_length(output) == contig_to_coding_length


def test_write_contig_coding_length_casts_float_to_int(tmp_path):
    output = tmp_path / "coding_length.tsv"
    io_manager.write_contig_coding_length({"contig1": 120.0}, output)
    assert output.read_text() == "contig1\t120\n"


def test_clear_previous_outputs(tmp_path):
    final_bins_dir = tmp_path / "final_bins"
    report_dir = tmp_path / "input_bins_quality_reports"
    final_bins_dir.mkdir()
    report_dir.mkdir()

    old_bin = final_bins_dir / "old_prefix_bin1.fa"
    old_report = report_dir / "input_bins_1.A.tsv"
    unrelated_file = final_bins_dir / "notes.txt"
    for path in (old_bin, old_report, unrelated_file):
        path.write_text("x")

    io_manager.clear_previous_outputs(final_bins_dir, report_dir)

    assert not old_bin.exists()
    assert not old_report.exists()
    assert unrelated_file.exists()


def test_clear_previous_outputs_missing_dirs(tmp_path):
    # Nothing to clear on a fresh output directory
    io_manager.clear_previous_outputs(
        tmp_path / "final_bins", tmp_path / "input_bins_quality_reports"
    )


def test_write_original_bin_metrics(bin1, bin2, tmp_path):
    # Test that `write_original_bin_metrics` correctly writes bin metrics to files

    temp_directory = tmp_path / "test_output"
    # Call the function with mock data
    io_manager.write_original_bin_metrics([bin1, bin2], temp_directory)

    # Check if the output directory was created
    assert temp_directory.exists(), "Output directory should be created."

    # Check that the correct files are created
    expected_files = [
        temp_directory / "input_bins_1.test1.tsv",
        temp_directory / "input_bins_2.test2.tsv",
    ]

    assert temp_directory.exists(), (
        f"Expected temp_directory {temp_directory} was not created."
    )

    for file in expected_files:
        assert file.exists(), f"Expected file {file} was not created."


def test_read_contig2bin_table(tmp_path):
    table = tmp_path / "final_contig_to_bin.tsv"
    table.write_text("contig1\tbin1\ncontig2\tbin1\ncontig3\tbin2\n")
    assert io_manager.read_contig2bin_table(table) == {
        "contig1": "bin1",
        "contig2": "bin1",
        "contig3": "bin2",
    }


def test_write_bins_fasta_from_contig_to_bin(tmp_path):
    contigs_fasta = tmp_path / "contigs.fasta"
    contigs_fasta.write_text(">contig1\nAAAA\n>contig2\nCCCC\n>contig3\nGGGG\n")
    outdir = tmp_path / "final_bins"
    outdir.mkdir()
    (outdir / "bin1.fa").write_text(">stale\nTTTT\n")

    io_manager.write_bins_fasta_from_contig_to_bin(
        {"contig1": "bin1", "contig3": "bin1", "contig2": "bin2"},
        contigs_fasta,
        outdir,
    )

    assert (outdir / "bin1.fa").read_text() == ">contig1\nAAAA\n>contig3\nGGGG\n"
    assert (outdir / "bin2.fa").read_text() == ">contig2\nCCCC\n"
