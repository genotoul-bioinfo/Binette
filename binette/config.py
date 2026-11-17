"""
Module      : config
Description : Parameters of a Binette run and paths of its output files.
Copyright   : (c) Jean Mainguy
License     : GPL-3.0
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from binette.checkpoint import CHECKPOINT_FILENAME


@dataclass(frozen=True)
class OutputPaths:
    """Paths of the files and directories written by Binette in its output directory."""

    outdir: Path
    tmp_dir: Path
    manifest: Path
    faa_file: Path
    coding_length_file: Path
    diamond_result_file: Path
    final_bin_report: Path
    final_contig2bin_table: Path
    original_bin_report_dir: Path
    final_bins_dir: Path
    index_to_contig_file: Path
    all_bin_compo_file: Path

    @classmethod
    def from_outdir(cls, outdir: Path) -> "OutputPaths":
        tmp_dir = outdir / "temporary_files"
        return cls(
            outdir=outdir,
            tmp_dir=tmp_dir,
            manifest=tmp_dir / CHECKPOINT_FILENAME,
            faa_file=tmp_dir / "assembly_proteins.faa.gz",
            coding_length_file=tmp_dir / "contig_coding_length.tsv.gz",
            diamond_result_file=tmp_dir / "diamond_result.tsv.gz",
            final_bin_report=outdir / "final_bins_quality_reports.tsv",
            final_contig2bin_table=outdir / "final_contig_to_bin.tsv",
            original_bin_report_dir=outdir / "input_bins_quality_reports",
            final_bins_dir=outdir / "final_bins",
            index_to_contig_file=outdir / "index_to_contig.tsv",
            all_bin_compo_file=outdir / "all_bins_quality_reports.tsv",
        )


@dataclass(frozen=True)
class RunConfig:
    """Parameters of a Binette run."""

    contigs: Path
    bin_dirs: list[Path] | None = None
    contig2bin_tables: list[Path] | None = None
    proteins: Path | None = None
    outdir: Path = Path("results")
    prefix: str = "binette"
    threads: int = 1
    min_completeness: int = 40
    max_contamination: int = 10
    min_length: int = 200_000
    max_length: int = 10_000_000
    contamination_weight: float = 2.0
    fasta_extensions: tuple[str, ...] = (".fasta", ".fa", ".fna")
    checkm2_db: Path | None = None
    low_mem: bool = False
    resume: bool = True
    force: bool = False
    debug: bool = False
    write_fasta_bins: bool = True
    show_progress: bool = True

    def selection_parameters(self) -> dict[str, Any]:
        """
        Parameters that change the final results without changing the intermediate
        results (proteins, DIAMOND hits). Writing bin FASTA files is not one of them:
        it is tracked separately.
        """
        return {
            "min_completeness": self.min_completeness,
            "max_contamination": self.max_contamination,
            "min_length": self.min_length,
            "max_length": self.max_length,
            "contamination_weight": self.contamination_weight,
            "prefix": self.prefix,
            "debug": self.debug,
        }
