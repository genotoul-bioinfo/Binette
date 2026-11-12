from collections import Counter
from pathlib import Path
from unittest.mock import patch

import pytest

from binette.config import OutputPaths
from binette.pipeline import (
    clean_temporary_files,
    get_diamond_hits,
    get_proteins,
    parse_input_files,
)

FAA_CONTENT = ">contig1_1\nMLKPACGT\n>contig2_1\nMMMKPTGCA\n>contig2_2\nMMMAAAA\n>contig3_1\nMLPALP\n"


CONTIGS_IN_BINS = {"contig1", "contig2", "contig3"}


@pytest.fixture
def contigs_fasta(tmp_path):
    fasta = tmp_path / "contigs.fasta"
    fasta.write_text(">contig1\nACGT\n>contig2\nACGT\n>contig3\nACGT\n")
    return fasta


def test_get_proteins_reuse_predicted(tmp_path, contigs_fasta):
    faa_file = tmp_path / "proteins.faa"
    faa_file.write_text(FAA_CONTENT)
    coding_length_file = tmp_path / "coding_length.tsv"
    coding_length_file.write_text("contig1\t3\ncontig2\t4\ncontig3\t2\n")

    with patch("binette.cds.predict") as mock_predict:
        contig_to_genes, contig_to_coding_len, outputs = get_proteins(
            faa_file=faa_file,
            coding_length_file=coding_length_file,
            contigs_fasta=contigs_fasta,
            contigs_in_bins=CONTIGS_IN_BINS,
            proteins=None,
            threads=1,
            reuse=True,
        )
        mock_predict.assert_not_called()

    assert len(contig_to_genes) == 3
    assert contig_to_genes["contig2"] == ["MMMKPTGCA", "MMMAAAA"]
    assert contig_to_coding_len == {"contig1": 3, "contig2": 4, "contig3": 2}
    assert outputs == [faa_file, coding_length_file]


def test_get_proteins_reuse_user_proteins(tmp_path, contigs_fasta):
    faa_file = tmp_path / "proteins.faa"
    faa_file.write_text(FAA_CONTENT)

    with patch("binette.cds.filter_faa_file") as mock_filter:
        contig_to_genes, contig_to_coding_len, outputs = get_proteins(
            faa_file=faa_file,
            coding_length_file=tmp_path / "coding_length.tsv",
            contigs_fasta=contigs_fasta,
            contigs_in_bins=CONTIGS_IN_BINS,
            proteins=tmp_path / "user_proteins.faa",
            threads=1,
            reuse=True,
        )
        mock_filter.assert_not_called()

    assert len(contig_to_genes) == 3
    assert contig_to_coding_len is None
    assert outputs == [faa_file]


def test_get_proteins_from_user_proteins(tmp_path, contigs_fasta):
    user_proteins = tmp_path / "user_proteins.faa"
    user_proteins.write_text(FAA_CONTENT + ">contig_not_in_bins_1\nMMMM\n")
    faa_file = tmp_path / "filtered_proteins.faa"

    contig_to_genes, contig_to_coding_len, outputs = get_proteins(
        faa_file=faa_file,
        coding_length_file=tmp_path / "coding_length.tsv",
        contigs_fasta=contigs_fasta,
        contigs_in_bins=CONTIGS_IN_BINS,
        proteins=user_proteins,
        threads=1,
        reuse=False,
    )

    assert set(contig_to_genes) == CONTIGS_IN_BINS
    assert contig_to_coding_len is None
    assert outputs == [faa_file]
    assert "contig_not_in_bins" not in faa_file.read_text()


def test_get_proteins_predict(tmp_path, contigs_fasta):
    faa_file = tmp_path / "proteins.faa.gz"
    coding_length_file = tmp_path / "coding_length.tsv"

    with patch("binette.cds.predict") as mock_predict:
        mock_predict.return_value = ({"contig1": ["MLK"]}, {"contig1": 9.0})
        contig_to_genes, contig_to_coding_len, outputs = get_proteins(
            faa_file=faa_file,
            coding_length_file=coding_length_file,
            contigs_fasta=contigs_fasta,
            contigs_in_bins={"contig1"},
            proteins=None,
            threads=2,
            reuse=False,
        )
        mock_predict.assert_called_once()

    assert contig_to_genes == {"contig1": ["MLK"]}
    assert contig_to_coding_len == {"contig1": 9.0}
    assert outputs == [faa_file, coding_length_file]
    assert coding_length_file.read_text() == "contig1\t9\n"


@pytest.mark.parametrize("reuse", [True, False])
def test_get_diamond_hits(tmp_path, contigs_fasta, reuse):
    faa_file = tmp_path / "proteins.faa"
    diamond_result_file = tmp_path / "diamond_result.tsv.gz"
    db = tmp_path / "db.dmnd"

    contig_to_kegg_id = {
        "contig1": Counter({"K12345": 1, "K67890": 1}),
        "contig2": Counter({"K23456": 1}),
    }

    with (
        patch(
            "binette.diamond.get_contig_to_kegg_id", return_value=contig_to_kegg_id
        ) as mock_parse,
        patch("binette.diamond.run") as mock_run,
    ):
        result = get_diamond_hits(
            faa_file=faa_file,
            contigs_fasta=contigs_fasta,
            contigs_in_bins=CONTIGS_IN_BINS,
            diamond_result_file=diamond_result_file,
            diamond_db_path=db,
            threads=3,
            low_mem=True,
            reuse=reuse,
        )

    assert result == contig_to_kegg_id
    mock_parse.assert_called_once_with(diamond_result_file.as_posix())
    if reuse:
        mock_run.assert_not_called()
    else:
        mock_run.assert_called_once_with(
            faa_file.as_posix(),
            diamond_result_file.as_posix(),
            db.as_posix(),
            (tmp_path / "diamond_result.log").as_posix(),
            3,
            low_mem=True,
        )


def test_get_diamond_hits_inconsistent_contigs(tmp_path, contigs_fasta):
    with (
        patch(
            "binette.diamond.get_contig_to_kegg_id",
            return_value={"unknown_contig": Counter({"K1": 1})},
        ),
        pytest.raises(AssertionError),
    ):
        get_diamond_hits(
            faa_file=tmp_path / "proteins.faa",
            contigs_fasta=contigs_fasta,
            contigs_in_bins=CONTIGS_IN_BINS,
            diamond_result_file=tmp_path / "diamond_result.tsv.gz",
            diamond_db_path=tmp_path / "db.dmnd",
            threads=1,
            low_mem=False,
            reuse=True,
        )


def test_parse_input_files_with_contig2bin_tables(tmp_path):
    bin_set1 = tmp_path / "bin_set1.tsv"
    bin_set1.write_text("contig1\tbin1A\ncontig2\tbin1B\n")
    bin_set2 = tmp_path / "bin_set2.tsv"
    bin_set2.write_text("contig3\tbin2A\ncontig4\ttbin2B\n")

    fasta_file = tmp_path / "assembly.fasta"
    fasta_file_content = ">contig1\nACGT\n>contig2\nTGCA\n>contig3\nAAAA\n>contig4\nCCCC\n>contig5\nCGTCGCT\n"
    fasta_file.write_text(fasta_file_content)

    # Call the function and capture the return values
    inputs = parse_input_files(None, [bin_set1, bin_set2], fasta_file)
    contig_key_to_bin = inputs.contig_key_to_original_bin
    contigs_in_bins = inputs.contigs_in_bins
    contig_id_to_length = inputs.contig_to_length

    # # Perform assertions on the returned values
    assert isinstance(contig_key_to_bin, dict)
    assert isinstance(contigs_in_bins, list)
    assert isinstance(contig_id_to_length, dict)

    assert len(contig_key_to_bin) == 4
    assert set(contigs_in_bins) == {"contig1", "contig2", "contig3", "contig4"}
    assert len(contig_id_to_length) == 4
    assert set(inputs.fingerprints) == {"contigs", "bins"}


def test_parse_input_files_fingerprints(tmp_path):
    fasta_file = tmp_path / "assembly.fasta"
    fasta_file.write_text(">contig1\nACGT\n>contig2\nTGCA\n>contig3\nAAAA\n")
    bin_table = tmp_path / "bin_set1.tsv"

    def fingerprints(table_content, fasta_content=None):
        bin_table.write_text(table_content)
        if fasta_content is not None:
            fasta_file.write_text(fasta_content)
        return parse_input_files(None, [bin_table], fasta_file).fingerprints

    reference = fingerprints("contig1\tbinA\ncontig2\tbinB\n")

    # Same input gives the same fingerprints
    assert fingerprints("contig1\tbinA\ncontig2\tbinB\n") == reference

    # Same contigs, different bin membership: only bins fingerprint changes
    moved = fingerprints("contig1\tbinA\ncontig2\tbinA\n")
    assert moved["contigs"] == reference["contigs"]
    assert moved["bins"] != reference["bins"]

    # Different set of contigs in bins
    assert (
        fingerprints("contig1\tbinA\ncontig3\tbinB\n")["contigs"]
        != (reference["contigs"])
    )

    # Same names, modified sequence
    modified = fingerprints(
        "contig1\tbinA\ncontig2\tbinB\n", ">contig1\nACGT\n>contig2\nTGCC\n"
    )
    assert modified["contigs"] != reference["contigs"]

    # Contig order in the FASTA file does not matter
    reordered = fingerprints(
        "contig1\tbinA\ncontig2\tbinB\n", ">contig2\nTGCA\n>contig1\nACGT\n"
    )
    assert reordered == reference


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
    inputs = parse_input_files(bin_dirs, contig2bin_tables, fasta_file)
    contig_key_to_bin = inputs.contig_key_to_original_bin
    contigs_in_bins = inputs.contigs_in_bins
    contig_id_to_length = inputs.contig_to_length

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


def test_clean_temporary_files(tmp_path):
    paths = OutputPaths.from_outdir(tmp_path / "results")
    paths.tmp_dir.mkdir(parents=True)
    for path in (
        paths.manifest,
        paths.faa_file,
        paths.coding_length_file,
        paths.diamond_result_file,
        paths.tmp_dir / "diamond_result.log",
    ):
        path.write_text("content")
    paths.final_bin_report.write_text("final results")

    clean_temporary_files(paths)

    # Only the checkpoint is kept in the temporary directory
    assert list(paths.tmp_dir.iterdir()) == [paths.manifest]
    # Final results are untouched
    assert paths.final_bin_report.read_text() == "final results"

    # Nothing left to clean: no error
    clean_temporary_files(paths)
