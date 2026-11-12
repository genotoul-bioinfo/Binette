#!/usr/bin/env python
"""
Module      : Main
Description : Command line interface of Binette.
Copyright   : (c) Jean Mainguy, 28 nov. 2022
License     : GPL-3.0
Maintainer  : Jean Mainguy
Portability : POSIX
"""

import logging
import sys
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.logging import RichHandler

import binette as binette_init
from binette import checkpoint, pipeline
from binette.config import RunConfig

logger = logging.getLogger(__name__)
err_console = Console(stderr=True)


def version_callback(
    value: bool,
    ctx: typer.Context,
):
    """Prints the version and exits if --version is passed."""
    if ctx.resilient_parsing:
        return

    if value:
        typer.echo(f"Binette {binette_init.__version__}")
        raise typer.Exit()


def setup_logging(verbose: bool = False, quiet: bool = False):
    """Sets up logging configuration based on verbosity flags."""
    if quiet and verbose:
        raise typer.BadParameter("Cannot specify both --verbose and --quiet")

    if quiet:
        lvl = logging.WARNING
    elif verbose:
        lvl = logging.DEBUG
    else:
        lvl = logging.INFO

    # Set up logging
    logging.basicConfig(
        level=lvl,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(console=err_console)],
    )

    # Only log startup messages if not in quiet mode
    if not quiet:
        logger.info("Program started")
        logger.info(f"Command line: {' '.join(sys.argv)}")


def verbose_callback(
    verbose: bool,
):
    """Sets the logging level to DEBUG if --verbose is passed."""
    # This is a placeholder - actual setup happens in the main function
    return verbose


def quiet_callback(
    quiet: bool,
):
    """Sets the logging level to WARNING if --quiet is passed."""
    # This is a placeholder - actual setup happens in the main function
    return quiet


def preprocess_args():
    """
    Typer doesn't support whitespace-separated multi-value options.

    We preprocess the sysargv so that:
    - python3 app.py some_command --filters filter1 filter2 filter3 \
      --environments env1 env2 env3

    becomes:
    - python3 app.py some_command --filters filter1 --filters filter2 --filters filter3 --environments env1 --environments env2 --environments env3

    """

    logger.debug(f"Initial CLI command is: {sys.argv}")

    # get main cmd
    final_cmd = []
    for _, arg in enumerate(sys.argv):
        if any(arg.startswith(_) for _ in ["-", "--"]):
            break
        else:
            final_cmd.append(arg)
    logger.debug(f"Main command is: {final_cmd}")

    # get options and their values
    for idx, arg in enumerate(sys.argv):
        if any(arg.startswith(_) for _ in ["-", "--"]):
            opt_values = []
            for value in sys.argv[idx + 1 :]:
                if any(value.startswith(_) for _ in ["-", "--"]):
                    break
                else:
                    opt_values.append(value)

            if len(opt_values) >= 1:
                [final_cmd.extend([arg, opt_value]) for opt_value in opt_values]
            else:
                final_cmd.append(arg)

    # replace by reformatted
    logger.debug(f"Final command is: {final_cmd}")
    sys.argv = final_cmd


def validate_bin_inputs(
    bin_dirs: list[Path] | None, contig2bin_tables: list[Path] | None
) -> None:
    """Validate that exactly one of bin_dirs or contig2bin_tables is provided."""
    if bin_dirs is None and contig2bin_tables is None:
        raise typer.BadParameter(
            "Error: Either --bin-dirs or --contig2bin_tables must be provided. None were given."
        )

    if bin_dirs is not None and contig2bin_tables is not None:
        raise typer.BadParameter(
            "Error: Either --bin-dirs or --contig2bin_tables must be provided, but not both."
        )


# Create the Typer app with no args help enabled and rich formatting
app = typer.Typer(
    name="binette",
    help=f"Binette: binning refinement tool to constructs high quality MAGs. Version: {binette_init.__version__}",
    add_completion=False,
    context_settings={"help_option_names": ["-h", "--help"]},
    rich_markup_mode="rich",
)


@app.command(
    help=f"Binette {binette_init.__version__}: fast and accurate binning refinement tool to constructs high quality MAGs from the output of multiple binning tools.",
    no_args_is_help=True,
)
def binette(
    # Input arguments - Mutually exclusive group (handled in code)
    bin_dirs: Annotated[
        list[Path] | None,
        typer.Option(
            "--bin_dirs",
            "-d",
            help="List of bin folders containing each bin in a fasta file.",
            # callback=lambda x: [is_valid_file(str(p)) for p in x] if x else None,
            exists=True,
            rich_help_panel="Input Arguments",
        ),
    ] = None,
    contig2bin_tables: Annotated[
        list[Path] | None,
        typer.Option(
            "--contig2bin_tables",
            "-b",
            help="List of contig2bin tables with two columns: contig, bin.",
            exists=True,
            rich_help_panel="Input Arguments",
        ),
    ] = None,
    contigs: Annotated[
        Path,
        typer.Option(
            "--contigs",
            "-c",
            help="Contigs in FASTA format.",
            exists=True,
            rich_help_panel="Input Arguments",
        ),
    ] = ...,  # Required
    proteins: Annotated[
        Path | None,
        typer.Option(
            "--proteins",
            "-p",
            help="FASTA file of predicted proteins in Prodigal format (>contigID_geneID). Skips the gene prediction step if provided.",
            exists=True,
            rich_help_panel="Input Arguments",
        ),
    ] = None,
    # Output & runtime control
    outdir: Annotated[
        Path,
        typer.Option(
            "--outdir",
            "-o",
            help="Output directory.",
            rich_help_panel="Output and Runtime Control",
        ),
    ] = Path("results"),
    prefix: Annotated[
        str,
        typer.Option(
            "--prefix",
            help="Prefix to add to final bin names (e.g. '--prefix sample1' will produce 'sample1_bin1.fa', 'sample1_bin2.fa').",
            rich_help_panel="Output and Runtime Control",
        ),
    ] = "binette",
    threads: Annotated[
        int,
        typer.Option(
            "--threads",
            "-t",
            help="Number of threads to use.",
            rich_help_panel="Output and Runtime Control",
        ),
    ] = 1,
    verbose: Annotated[
        bool,
        typer.Option(
            "--verbose",
            "-v",
            help="Enable verbose mode (show detailed debug information).",
            callback=verbose_callback,
            rich_help_panel="Output and Runtime Control",
        ),
    ] = False,
    quiet: Annotated[
        bool,
        typer.Option(
            "--quiet",
            "-q",
            help="Enable quiet mode (only show warnings and errors).",
            callback=quiet_callback,
            rich_help_panel="Output and Runtime Control",
        ),
    ] = False,
    debug: Annotated[
        bool,
        typer.Option(
            help="Activate debug mode.",
            rich_help_panel="Output and Runtime Control",
        ),
    ] = False,
    # Bin filtering & scoring
    min_completeness: Annotated[
        int,
        typer.Option(
            "--min_completeness",
            help="Minimum completeness required for intermediate bin creation and final bin selection.",
            rich_help_panel="Bin Filtering and Scoring",
        ),
    ] = 40,
    max_contamination: Annotated[
        int,
        typer.Option(
            "--max_contamination",
            help="Maximum contamination allowed for intermediate bin creation and final bin selection.",
            rich_help_panel="Bin Filtering and Scoring",
        ),
    ] = 10,
    min_length: Annotated[
        int,
        typer.Option(
            "--min_length",
            help="Minimum length (bp) required for intermediate bin creation and final bin selection.",
            rich_help_panel="Bin Filtering and Scoring",
        ),
    ] = 200_000,
    max_length: Annotated[
        int,
        typer.Option(
            "--max_length",
            help="Maximum length (bp) allowed for intermediate bin creation and final bin selection.",
            rich_help_panel="Bin Filtering and Scoring",
        ),
    ] = 10_000_000,
    contamination_weight: Annotated[
        float,
        typer.Option(
            "--contamination_weight",
            "-w",
            help="Bins are scored as: completeness - weight * contamination. A lower weight favors completeness over low contamination.",
            rich_help_panel="Bin Filtering and Scoring",
        ),
    ] = 2.0,
    # Advanced options
    fasta_extensions: Annotated[
        list[str],
        typer.Option(
            "--fasta_extensions",
            "-e",
            help="FASTA file extensions to search for in bin directories (used with --bin_dirs).",
            rich_help_panel="Advanced Options",
        ),
    ] = [  # noqa: B006
        ".fasta",
        ".fa",
        ".fna",
    ],
    checkm2_db: Annotated[
        Path | None,
        typer.Option(
            "--checkm2_db",
            help="Path to CheckM2 diamond database. By default the database set via <checkm2 database> is used.",
            rich_help_panel="Advanced Options",
        ),
    ] = None,
    low_mem: Annotated[
        bool,
        typer.Option(
            "--low_mem",
            help="Enable low-memory mode for Diamond.",
            rich_help_panel="Advanced Options",
        ),
    ] = False,
    resume: Annotated[
        bool,
        typer.Option(
            help=(
                "Reuse intermediate results (proteins, DIAMOND hits) "
                "from a previous run with the same inputs "
                "(disable with --no-resume)."
            ),
            rich_help_panel="Advanced Options",
        ),
    ] = True,
    force: Annotated[
        bool,
        typer.Option(
            "--force",
            help=(
                "Allow using an output directory that contains results of a run with "
                "different inputs. Previous intermediate "
                "files are discarded and results are overwritten."
            ),
            rich_help_panel="Advanced Options",
        ),
    ] = False,
    version: Annotated[
        bool,
        typer.Option(
            "--version",
            help="Show version and exit.",
            callback=version_callback,
        ),
    ] = None,
    progress: Annotated[
        bool,
        typer.Option(
            help="Show progress bar while fetching pangenomes (disable with --no-progress).",
            rich_help_panel="Output and Runtime Control",
        ),
    ] = True,
    write_fasta_bins: Annotated[
        bool,
        typer.Option(
            help="Write final selected bins as FASTA files (disable with --no-write-fasta-bins).",
            rich_help_panel="Output and Runtime Control",
        ),
    ] = True,
) -> int:
    """Orchestrate the execution of the program"""

    # Set up logging based on verbosity flags
    setup_logging(verbose=verbose, quiet=quiet)

    validate_bin_inputs(bin_dirs, contig2bin_tables)

    config = RunConfig(
        contigs=contigs,
        bin_dirs=bin_dirs,
        contig2bin_tables=contig2bin_tables,
        proteins=proteins,
        outdir=outdir,
        prefix=prefix,
        threads=threads,
        min_completeness=min_completeness,
        max_contamination=max_contamination,
        min_length=min_length,
        max_length=max_length,
        contamination_weight=contamination_weight,
        fasta_extensions=tuple(fasta_extensions),
        checkm2_db=checkm2_db,
        low_mem=low_mem,
        resume=resume,
        force=force,
        debug=debug,
        write_fasta_bins=write_fasta_bins,
        show_progress=progress and not quiet,
    )

    try:
        return pipeline.run(config)
    except checkpoint.OutputDirConflictError as exc:
        logger.error(str(exc))
        raise typer.Exit(code=1) from exc


def main():
    preprocess_args()

    app()
