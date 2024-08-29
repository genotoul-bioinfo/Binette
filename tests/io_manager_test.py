import pytest
from binette import io_manager
from pathlib import Path
from unittest.mock import patch




class Bin:
    def __init__(self, bin_id, origin, name, completeness, contamination, score, length, N50, contigs):
        self.id = bin_id
        self.origin = {origin}
        self.name = name
        self.completeness = completeness
        self.contamination = contamination
        self.score = score
        self.length = length
        self.N50 = N50
        self.contigs = contigs

@pytest.fixture
def bin1():
    return Bin(1, 'origin1', 'name1', 90, 5, 80, 1000, 500, ['contig1', 'contig3'])

@pytest.fixture
def bin2():
    return Bin(2, 'origin2', 'name2', 85, 8, 75, 1200, 600, ['contig2', 'contig4'])


def test_infer_bin_name_from_bin_inputs():
    # Mock input data
    input_bins = [
        '/path/to/bin1',
        '/path/to/bin2',
        '/path/to/bin3'
    ]

    # Call the function
    result = io_manager.infer_bin_set_names_from_input_paths(list(map(Path, input_bins)))

    # Define the expected output
    expected_result = {
        'bin1': Path('/path/to/bin1'),
        'bin2': Path('/path/to/bin2'),
        'bin3': Path('/path/to/bin3')
    }

    # Check if the output matches the expected dictionary
    assert result == expected_result

def test_infer_bin_name_from_single_path():
    # Mock input data
    input_bins = [
        '/path/to/bin1',
    ]

    # Call the function
    result = io_manager.infer_bin_set_names_from_input_paths(list(map(Path, input_bins)))

    # Define the expected output
    expected_result = {
        '/path/to/bin1': Path('/path/to/bin1'),
    }

    # Check if the output matches the expected dictionary
    assert result == expected_result


def test_infer_bin_name_from_bin_table_inputs():
    # Mock input data
    input_bins = [
        '/path/to/bin1.tsv',
        '/path/to/bin2.tsv',
        '/path/to/bin3.tsv'
    ]

    # Call the function
    result = io_manager.infer_bin_set_names_from_input_paths(list(map(Path, input_bins)))

    # Define the expected output
    expected_result = {
        'bin1': Path('/path/to/bin1.tsv'),
        'bin2': Path('/path/to/bin2.tsv'),
        'bin3': Path('/path/to/bin3.tsv')
    }

    # Check if the output matches the expected dictionary
    assert result == expected_result


def test_infer_bin_name_from_bin_table_with_different_ext():
    # Mock input data
    input_bins = [
        '/path/to/bin1.tsv',
        '/path/to/bin2.tsv',
        '/path/to/bin3.txt'
    ]

    # Call the function
    result = io_manager.infer_bin_set_names_from_input_paths(list(map(Path, input_bins)))

    # Define the expected output
    expected_result = {
        'bin1.tsv': Path('/path/to/bin1.tsv'),
        'bin2.tsv': Path('/path/to/bin2.tsv'),
        'bin3.txt': Path('/path/to/bin3.txt')
    }

    # Check if the output matches the expected dictionary
    assert result == expected_result

def test_infer_bin_name_from_bin_table_with_different_dir():
    # Mock input data
    input_bins = [
        '/path/to/bins',
        '/path2/result_bins',
        '/path2/result/bins',
    ]

    # Call the function
    result = io_manager.infer_bin_set_names_from_input_paths(list(map(Path, input_bins)))

    # Define the expected output
    expected_result = {
        'path/to/bins' : Path('/path/to/bins'),
        'path2/result_bins': Path('/path2/result_bins'),
        'path2/result/bins': Path('/path2/result/bins'),
    }

    # Check if the output matches the expected dictionary
    assert result == expected_result
    
def test_get_paths_common_prefix_suffix():
    # Test case 1: No paths provided
    assert io_manager.get_paths_common_prefix_suffix([]) == ([], [], [])

    # # Test case 2: Single path
    assert io_manager.get_paths_common_prefix_suffix([Path('/home/user/project')]) == (['/', 'home', 'user', 'project'], ['/', 'home', 'user', 'project'], [])

    # Test case 3: Multiple paths with common prefix and suffix
    paths = [Path('/home/user/project/src'), Path('/home/user/project/docs'), Path('/home/user/project/tests')]
    assert io_manager.get_paths_common_prefix_suffix(paths) == (['/', 'home', 'user', 'project'], [], [])

    # Test case 4: Multiple paths with no common prefix or suffix
    paths = [Path('/var/log/syslog'), Path('/usr/local/bin/python'), Path('/etc/nginx/nginx.conf')]
    assert io_manager.get_paths_common_prefix_suffix(paths) == (['/'], [], [])

    # Test case 5: Multiple paths with common suffix
    paths = [Path('/home/user/docs/report.txt'), Path('/home/admin/docs/report.txt')]
    assert io_manager.get_paths_common_prefix_suffix(paths) == (['/', 'home'], ['docs', 'report.txt'], ['.txt'])

    # Test case 6: Paths with a deeper common prefix and suffix
    paths = [Path('/data/project_a/results/output.txt'), Path('/data/project_b/results/output.txt')]
    assert io_manager.get_paths_common_prefix_suffix(paths) == (['/', 'data'], ['results', 'output.txt'], ['.txt'])

    # Test case 7: Paths with only the root as common prefix and different suffix
    paths = [Path('/project_a/output.txt'), Path('/project_b/output.txt')]
    assert io_manager.get_paths_common_prefix_suffix(paths) == (['/'], ['output.txt'], ['.txt'])

    # Test case 8: Paths with only the root as common prefix and different suffix
    paths = [Path('/project_a/output.txt'), Path('/project_a/output.tsv')]
    assert io_manager.get_paths_common_prefix_suffix(paths) == (['/', 'project_a'], [], [])

def test_write_bin_info(tmp_path, bin1, bin2):
    # Mock input data
    bins = [bin1, bin2]

    output_file = tmp_path / "output.tsv"

    # Call the function
    io_manager.write_bin_info(bins, output_file)

    # Check if the file was created and its content matches the expected output
    assert Path(output_file).exists()

    with open(output_file, "r") as f:
        content = f.read()
        assert "bin_id\torigin\tname\tcompleteness\tcontamination\tscore\tsize\tN50\tcontig_count" in content
        assert "1\torigin1\tname1\t90\t5\t80\t1000\t500\t2" in content
        assert "2\torigin2\tname2\t85\t8\t75\t1200\t600\t2" in content


def test_write_bin_info_add_contig(tmp_path, bin1, bin2):
    # Mock input data
    bins = [bin1, bin2]

    output_file = tmp_path / "output.tsv"

    # Call the function
    io_manager.write_bin_info(bins, output_file, add_contigs=True)

    # Check if the file was created and its content matches the expected output
    assert Path(output_file).exists()

    with open(output_file, "r") as f:
        content = f.read()
        assert "bin_id\torigin\tname\tcompleteness\tcontamination\tscore\tsize\tN50\tcontig_count\tcontigs" in content
        assert "1\torigin1\tname1\t90\t5\t80\t1000\t500\t2\tcontig1;contig3" in content
        assert "2\torigin2\tname2\t85\t8\t75\t1200\t600\t2\tcontig2;contig4" in content




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
    io_manager.write_bins_fasta(selected_bins, contigs_fasta, outdir)

    # Check if the files were created and their content matches the expected output
    assert (outdir / "bin_1.fa").exists()
    assert (outdir / "bin_2.fa").exists()

    with open(outdir / "bin_1.fa", "r") as bin1_file:
        assert bin1_file.read() == ">contig1\nACGT\n>contig3\nAAAA\n"

    with open(outdir / "bin_2.fa", "r") as bin2_file:
        assert bin2_file.read() == ">contig2\nTGCA\n>contig4\nCCCC\n"
    

def test_check_contig_consistency_error():
    # Mock input data
    contigs_from_assembly = ['contig1', 'contig2', 'contig3']
    contigs_from_bins = ['contig2', 'contig3', 'contig4']
    assembly_file = "assembly.fasta"
    elsewhere_file = "external.fasta"

    with pytest.raises(AssertionError):
        # Call the function
        io_manager.check_contig_consistency(
            contigs_from_assembly, contigs_from_bins, assembly_file, elsewhere_file
        )

def test_check_contig_consistency_no_error():
    # Mock input data
    contigs_from_assembly = ['contig1', 'contig2', 'contig3', 'contig4']
    contigs_from_bins = ['contig1', 'contig2', 'contig3']
    assembly_file = "assembly.fasta"
    elsewhere_file = "external.fasta"

    io_manager.check_contig_consistency(
            contigs_from_assembly, contigs_from_bins, assembly_file, elsewhere_file
        )

@pytest.fixture
def temp_files(tmp_path):
    # Create temporary files for testing
    faa_file = tmp_path / "test_protein.faa"
    diamond_result_file = tmp_path / "test_diamond_result.txt"
    faa_file.touch()
    diamond_result_file.touch()
    yield str(faa_file), str(diamond_result_file)

def test_check_resume_file_exists(temp_files, caplog):
    # Test when both files exist
    faa_file, diamond_result_file = temp_files
    io_manager.check_resume_file(Path(faa_file), Path(diamond_result_file))
    assert "Protein file" not in caplog.text
    assert "Diamond result file" not in caplog.text

def test_check_resume_file_missing_faa(temp_files, caplog):
    # Test when faa_file is missing
    _, diamond_result_file = temp_files
    with pytest.raises(FileNotFoundError):
        io_manager.check_resume_file(Path("nonexistent.faa"), Path(diamond_result_file))
    assert "Protein file" in caplog.text
    assert "Diamond result file" not in caplog.text

def test_check_resume_file_missing_diamond(temp_files, caplog):
    # Test when diamond_result_file is missing
    faa_file, _ = temp_files
    with pytest.raises(FileNotFoundError):
        io_manager.check_resume_file(Path(faa_file), Path("nonexistent_diamond_result.txt"))
    assert "Protein file" not in caplog.text
    assert "Diamond result file" in caplog.text


@patch('binette.io_manager.write_bin_info')
def test_write_original_bin_metrics(mock_write_bin_info, bin1,bin2, tmp_path):
    # Test that `write_original_bin_metrics` correctly writes bin metrics to files

    temp_directory =  tmp_path / "test_output"

    mock_bins = {"set1":{bin1},
                 "set2":{bin2}}
    # Call the function with mock data
    io_manager.write_original_bin_metrics(mock_bins, temp_directory)

    # Check if the output directory was created
    assert temp_directory.exists(), "Output directory should be created."

    # Check that the correct files are created
    expected_files = [
        temp_directory / "input_bins_1.set1.tsv",
        temp_directory / "input_bins_2.set2.tsv"
    ]

    assert temp_directory.exists(), f"Expected temp_directory {temp_directory} was not created."

    # Check if `write_bin_info` was called correctly
    assert mock_write_bin_info.call_count == 2, "write_bin_info should be called once for each bin set."

    # Verify the specific calls to `write_bin_info`
    mock_write_bin_info.assert_any_call(mock_bins['set1'], expected_files[0])
    mock_write_bin_info.assert_any_call(mock_bins['set2'], expected_files[1])