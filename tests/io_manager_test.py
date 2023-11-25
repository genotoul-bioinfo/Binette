import pytest
from binette import io_manager
from pathlib import Path





class Bin:
    def __init__(self, bin_id, origin, name, completeness, contamination, score, length, N50, contigs):
        self.id = bin_id
        self.origin = origin
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
    result = io_manager.infer_bin_name_from_bin_inputs(input_bins)

    # Define the expected output
    expected_result = {
        '1': '/path/to/bin1',
        '2': '/path/to/bin2',
        '3': '/path/to/bin3'
    }

    # Check if the output matches the expected dictionary
    assert result == expected_result


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
    io_manager.write_bins_fasta(selected_bins, str(contigs_fasta), str(outdir))

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
