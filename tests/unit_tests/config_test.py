from dataclasses import fields
from pathlib import Path

from binette.checkpoint import BINETTE_OUTPUT_MARKERS, CHECKPOINT_FILENAME
from binette.config import OutputPaths, RunConfig


def test_output_paths_from_outdir():
    outdir = Path("results")
    paths = OutputPaths.from_outdir(outdir)

    assert paths.manifest == outdir / "temporary_files" / CHECKPOINT_FILENAME
    for field in fields(OutputPaths):
        path = getattr(paths, field.name)
        assert path == outdir or outdir in path.parents

    # Output files recognised by the output directory guard are the ones written
    written = {
        getattr(paths, f.name).relative_to(outdir).as_posix() for f in fields(paths)
    }
    assert set(BINETTE_OUTPUT_MARKERS) <= written


def test_selection_parameters_exclude_runtime_options():
    base = RunConfig(contigs=Path("contigs.fa"))
    runtime_changes = RunConfig(
        contigs=Path("contigs.fa"),
        threads=8,
        low_mem=True,
        resume=False,
        force=True,
        write_fasta_bins=False,
        show_progress=False,
    )
    assert base.selection_parameters() == runtime_changes.selection_parameters()

    assert (
        RunConfig(contigs=Path("contigs.fa"), prefix="other").selection_parameters()
        != base.selection_parameters()
    )
