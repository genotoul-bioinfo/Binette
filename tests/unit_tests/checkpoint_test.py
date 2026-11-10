import json
import logging
from pathlib import Path

import pytest

from binette import checkpoint
from binette.checkpoint import Checkpoint, OutputDirConflictError, check_outdir

RUN_INPUTS = {
    "contigs": "c" * 64,
    "bins": "b" * 64,
    "proteins": None,
    "checkm2_db": {"name": "db.dmnd", "size": 10},
}


@pytest.fixture
def outdir(tmp_path) -> Path:
    return tmp_path / "results"


@pytest.fixture
def manifest_path(outdir) -> Path:
    return outdir / "temporary_files" / checkpoint.CHECKPOINT_FILENAME


def write(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


def test_sha256_file(tmp_path):
    path = write(tmp_path / "f.txt", "hello")
    assert (
        checkpoint.sha256_file(path)
        == "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
    )


def test_bins_fingerprint_order_independent():
    bins_a = {
        "A": [
            {"bin_name": "bin1", "contigs": ["c1", "c2"]},
            {"bin_name": "bin2", "contigs": ["c3"]},
        ]
    }
    bins_b = {
        "A": [
            {"bin_name": "bin2", "contigs": ["c3"]},
            {"bin_name": "bin1", "contigs": ["c2", "c1"]},
        ]
    }
    assert checkpoint.bins_fingerprint(bins_a) == checkpoint.bins_fingerprint(bins_b)

    bins_renamed = {"B": bins_a["A"]}
    assert checkpoint.bins_fingerprint(bins_a) != checkpoint.bins_fingerprint(
        bins_renamed
    )


def test_checkm2_db_identity_ignores_location(tmp_path):
    db1 = write(tmp_path / "a" / "db.dmnd", "database")
    db2 = write(tmp_path / "b" / "db.dmnd", "database")
    assert checkpoint.checkm2_db_identity(db1) == checkpoint.checkm2_db_identity(db2)


def test_record_and_reload(outdir, manifest_path):
    output = write(outdir / "temporary_files" / "proteins.faa", ">p\nMMM\n")
    ckpt = Checkpoint.new(manifest_path, outdir, RUN_INPUTS)
    ckpt.record_step("proteins", {"x": 1}, [output])

    assert manifest_path.exists()
    assert not checkpoint.partial_path(manifest_path).exists()

    reloaded, reason = Checkpoint.load(manifest_path, outdir)
    assert reason == ""
    assert reloaded.run_inputs == RUN_INPUTS
    assert reloaded.step_status("proteins", {"x": 1}) == (True, "")
    assert reloaded.step_output_digest("proteins", output) == checkpoint.sha256_file(
        output
    )

    # Output paths are stored relative to the root directory
    data = json.loads(manifest_path.read_text())
    assert list(data["steps"]["proteins"]["outputs"]) == [
        "temporary_files/proteins.faa"
    ]


def test_step_status_invalid_cases(outdir, manifest_path):
    output = write(outdir / "out.txt", "complete content")
    ckpt = Checkpoint.new(manifest_path, outdir, RUN_INPUTS)
    ckpt.record_step("step", {"x": 1}, [output])

    is_valid, reason = ckpt.step_status("unknown", {"x": 1})
    assert not is_valid and "no record" in reason

    is_valid, reason = ckpt.step_status("step", {"x": 2})
    assert not is_valid and "inputs changed" in reason

    # Truncated file (e.g. job killed while writing)
    output.write_text("complete")
    is_valid, reason = ckpt.step_status("step", {"x": 1})
    assert not is_valid and "incomplete" in reason

    # Same size, different content
    output.write_text("COMPLETE content")
    is_valid, reason = ckpt.step_status("step", {"x": 1})
    assert not is_valid and "changed" in reason

    output.unlink()
    is_valid, reason = ckpt.step_status("step", {"x": 1})
    assert not is_valid and "missing" in reason


def test_invalidate(outdir, manifest_path):
    output = write(outdir / "out.txt", "x")
    ckpt = Checkpoint.new(manifest_path, outdir, RUN_INPUTS)
    ckpt.record_step("a", {}, [output])
    ckpt.record_step("b", {}, [output])

    ckpt.invalidate("a", "not_recorded")

    reloaded, _ = Checkpoint.load(manifest_path, outdir)
    assert set(reloaded.data["steps"]) == {"b"}


@pytest.mark.parametrize(
    "content, expected_reason",
    [
        ("{not json", "unreadable"),
        (json.dumps({"format_version": 999}), "unsupported format version"),
        (json.dumps([1, 2]), "unsupported format version"),
    ],
)
def test_load_bad_manifest(outdir, manifest_path, content, expected_reason):
    write(manifest_path, content)
    ckpt, reason = Checkpoint.load(manifest_path, outdir)
    assert ckpt is None
    assert expected_reason in reason


def test_load_missing_manifest(outdir, manifest_path):
    ckpt, reason = Checkpoint.load(manifest_path, outdir)
    assert ckpt is None
    assert "no checkpoint file" in reason


def test_check_outdir_new(outdir, manifest_path):
    ckpt = check_outdir(outdir, manifest_path, RUN_INPUTS, force=False)
    assert ckpt.data["steps"] == {}
    assert manifest_path.exists()


def test_check_outdir_existing_unrelated_files(outdir, manifest_path):
    # An existing output directory without Binette results can be used
    write(outdir / "slurm.log", "log")
    ckpt = check_outdir(outdir, manifest_path, RUN_INPUTS, force=False)
    assert ckpt.run_inputs == RUN_INPUTS


def test_check_outdir_same_inputs(outdir, manifest_path):
    output = write(outdir / "out.txt", "x")
    first = check_outdir(outdir, manifest_path, RUN_INPUTS, force=False)
    first.record_step("step", {}, [output])

    second = check_outdir(outdir, manifest_path, dict(RUN_INPUTS), force=False)
    assert "step" in second.data["steps"]


def test_check_outdir_different_inputs(outdir, manifest_path, caplog):
    output = write(outdir / "out.txt", "x")
    check_outdir(outdir, manifest_path, RUN_INPUTS, force=False).record_step(
        "step", {}, [output]
    )
    other_inputs = RUN_INPUTS | {"bins": "z" * 64}

    with pytest.raises(OutputDirConflictError) as exc_info:
        check_outdir(outdir, manifest_path, other_inputs, force=False)
    assert "input bin sets" in str(exc_info.value)
    assert "contig sequences" not in str(exc_info.value)
    assert "--force" in str(exc_info.value)

    caplog.set_level(logging.WARNING)
    ckpt = check_outdir(outdir, manifest_path, other_inputs, force=True)
    assert ckpt.data["steps"] == {}
    assert ckpt.run_inputs == other_inputs
    assert "--force is set" in caplog.text

    # The directory now belongs to the new run
    reloaded, _ = Checkpoint.load(manifest_path, outdir)
    assert reloaded.run_inputs == other_inputs


@pytest.mark.parametrize(
    "legacy_output",
    [
        "final_bins_quality_reports.tsv",
        "temporary_files/diamond_result.tsv.gz",
        "final_bins/bin1.fa",
    ],
)
def test_check_outdir_legacy_results(outdir, manifest_path, legacy_output):
    write(outdir / legacy_output, "old")

    with pytest.raises(OutputDirConflictError) as exc_info:
        check_outdir(outdir, manifest_path, RUN_INPUTS, force=False)
    assert "cannot be verified" in str(exc_info.value)

    ckpt = check_outdir(outdir, manifest_path, RUN_INPUTS, force=True)
    assert ckpt.run_inputs == RUN_INPUTS


def test_check_outdir_corrupted_manifest(outdir, manifest_path):
    write(manifest_path, "{corrupted")

    with pytest.raises(OutputDirConflictError) as exc_info:
        check_outdir(outdir, manifest_path, RUN_INPUTS, force=False)
    assert "unreadable" in str(exc_info.value)

    ckpt = check_outdir(outdir, manifest_path, RUN_INPUTS, force=True)
    assert ckpt.run_inputs == RUN_INPUTS
