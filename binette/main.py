#!/usr/bin/env python
"""
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Jean Mainguy, 28 nov. 2022
License     : GPL-3.0
Maintainer  : Jean Mainguy
Portability : POSIX
"""

import logging
import os
import sys
from pathlib import Path
from typing import Annotated

import pyfastx
import typer
from rich.console import Console
from rich.logging import RichHandler

import binette as binette_init
from binette import bin_manager, bin_quality, cds, contig_manager, diamond
from binette import io_manager as io

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


# Create the Typer app with no args help enabled and rich formatting
app = typer.Typer(
    name="binette",
    help=f"Binette: binning refinement tool to constructs high quality MAGs. Version: {binette_init.__version__}",
    add_completion=False,
    context_settings={"help_option_names": ["-h", "--help"]},
    rich_markup_mode="rich",
)


def parse_input_files(
    bin_dirs: list[Path],
    contig2bin_tables: list[Path],
    contigs_fasta: Path,
    fasta_extensions: set[str] | None = None,
):
    """
    Parses input files to retrieve information related to bins and contigs.

    :param bin_dirs: List of paths to directories containing bin FASTA files.
    :param contig2bin_tables: List of paths to contig-to-bin tables.
    :param contigs_fasta: Path to the contigs FASTA file.
    :param temporary_dir: Path to the temporary directory to store intermediate files.
    :param fasta_extensions: Possible fasta extensions to look for in the bin directory.

    :return: A tuple containing:
        - List of original bins.
        - Dictionary mapping bins to lists of contigs.
        - Dictionary mapping contig names to their lengths.
    """
    if fasta_extensions is None:
        fasta_extensions = {".fasta", ".fa", ".fna"}

    if bin_dirs:
        logger.info("Parsing bin directories")
        bin_name_to_bin_dir = io.infer_bin_set_names_from_input_paths(bin_dirs)
        bin_set_name_to_bins_info = bin_manager.parse_bin_directories(
            bin_name_to_bin_dir, fasta_extensions
        )
    else:
        logger.info("Parsing bin2contig files")
        bin_name_to_bin_table = io.infer_bin_set_names_from_input_paths(
            contig2bin_tables
        )
        bin_set_name_to_bins_info = bin_manager.parse_contig2bin_tables(
            bin_name_to_bin_table
        )

    logger.info(f"Processing {len(bin_set_name_to_bins_info)} bin sets")
    for bin_set_id, bins_info in bin_set_name_to_bins_info.items():
        logger.info(f"  {bin_set_id} - {len(bins_info)} bins")

    contigs_in_bins = bin_manager.get_contigs_in_bin_sets(bin_set_name_to_bins_info)
    logger.info(f"Found {len(contigs_in_bins)} contigs in input bins")

    contig_to_index = contig_manager.make_contig_index(contigs_in_bins)

    contig_key_to_bin = bin_manager.make_bins_from_bins_info(
        bin_set_name_to_bins_info, contig_to_index, are_original_bins=True
    )

    # original_bins = bin_manager.dereplicate_bin_sets(bin_set_name_to_bins.values())

    logger.info(
        f"Parsing contig fasta file '{contigs_fasta}' to retrieve contig lengths"
    )

    contigs_in_bins_set = set(contigs_in_bins)
    contig_to_length = {
        name: len(seq)
        for name, seq in pyfastx.Fastx(contigs_fasta.as_posix())
        if name in contigs_in_bins_set
    }

    logger.debug("Finished parsing contig fasta file")
    # check if all contigs from input bins are present in contigs file
    unexpected_contigs = {
        contig for contig in contigs_in_bins if contig not in contig_to_length
    }

    if len(unexpected_contigs):
        raise ValueError(
            f"{len(unexpected_contigs)} contigs from the input bins were not "
            f"found in the contigs file '{contigs_fasta}'. "
            f"The missing contigs are: {', '.join(unexpected_contigs)}. "
            f"Please ensure all contigs from input bins are present in "
            f"contig file."
        )
    logger.debug("No unexpected contigs found")

    contig_id_to_length = {
        contig_to_index[name]: length for name, length in contig_to_length.items()
    }

    return (
        contig_key_to_bin,
        contigs_in_bins,
        contig_id_to_length,
        contig_to_index,
    )


def manage_protein_alignement(
    faa_file: Path,
    contigs_fasta: Path,
    contigs_in_bins: set[str],
    diamond_result_file: Path,
    checkm2_db: Path | None,
    threads: int,
    use_existing_protein_file: bool,
    resume_diamond: bool,
    low_mem: bool,
) -> tuple[dict[str, int], dict[str, list[str]], dict[str, int | None] | None]:
    """
    Predicts or reuses proteins prediction and runs diamond on them.

    :param faa_file: The path to the .faa file.
    :param contigs_fasta: The path to the contigs FASTA file.
    :param contigs_in_bins: Dictionary mapping bin names to lists of contigs.
    :param diamond_result_file: The path to the diamond result file.
    :param checkm2_db: The path to the CheckM2 database.
    :param threads: Number of threads for parallel processing.
    :param use_existing_protein_file: Boolean indicating whether to use an existing protein file.
    :param resume_diamond: Boolean indicating whether to resume diamond alignement.
    :param low_mem: Boolean indicating whether to use low memory mode.

    :return: A tuple containing dictionaries - contig_to_kegg_counter, contig_to_genes, and contig_to_coding_len.
    """

    # Predict or reuse proteins prediction and run diamond on them
    if use_existing_protein_file:
        logger.info(f"Parsing protein sequences from '{faa_file}'")
        contig_to_genes = cds.parse_faa_file(faa_file.as_posix())
        io.check_contig_consistency(
            contigs_in_bins,
            contig_to_genes,
            contigs_fasta.as_posix(),
            faa_file.as_posix(),
        )
        contig_to_coding_len = None
        logger.info(
            "Coding density will not be computed (using provided protein sequences)"
        )

    else:
        contigs_iterator = (
            (name, seq)
            for name, seq in pyfastx.Fastx(contigs_fasta.as_posix())
            if name in contigs_in_bins
        )
        contig_to_genes, contig_to_coding_len = cds.predict(
            contigs_iterator, faa_file.as_posix(), threads
        )
        logger.info("Coding density will be computed from freshly identified genes")

    if not resume_diamond:
        if checkm2_db is None:
            # get checkm2 db stored in checkm2 install
            diamond_db_path = diamond.get_checkm2_db()
        elif checkm2_db.exists():
            diamond_db_path = checkm2_db.as_posix()
        else:
            raise FileNotFoundError(checkm2_db)

        diamond_log = (
            diamond_result_file.parents[0]
            / f"{diamond_result_file.stem.split('.')[0]}.log"
        )

        diamond.run(
            faa_file.as_posix(),
            diamond_result_file.as_posix(),
            diamond_db_path,
            diamond_log.as_posix(),
            threads,
            low_mem=low_mem,
        )

    logger.info("Parsing diamond results")
    contig_to_kegg_counter = diamond.get_contig_to_kegg_id(
        diamond_result_file.as_posix()
    )

    # Check contigs from diamond vs input assembly consistency
    io.check_contig_consistency(
        contigs_in_bins,
        contig_to_kegg_counter,
        contigs_fasta.as_posix(),
        diamond_result_file.as_posix(),
    )

    return contig_to_kegg_counter, contig_to_genes, contig_to_coding_len


def log_selected_bin_info(
    selected_bins: list[bin_manager.Bin],
    hq_min_completeness: float,
    hq_max_conta: float,
):
    """
    Log information about selected bins based on quality thresholds.

    :param selected_bins: List of Bin objects to analyze.
    :param hq_min_completeness: Minimum completeness threshold for high-quality bins.
    :param hq_max_conta: Maximum contamination threshold for high-quality bins.

    This function logs information about selected bins that meet specified quality thresholds.
    It counts the number of high-quality bins based on completeness and contamination values.
    """

    # Log completeness and contamination in debug log
    logger.debug("High quality bins:")
    for sb in selected_bins:
        if sb.is_high_quality(
            min_completeness=hq_min_completeness, max_contamination=hq_max_conta
        ):
            logger.debug(
                f"  {sb} completeness={sb.completeness}, contamination={sb.contamination}"
            )

    # Count high-quality bins and single-contig high-quality bins
    hq_bins = len(
        [
            sb
            for sb in selected_bins
            if sb.is_high_quality(
                min_completeness=hq_min_completeness, max_contamination=hq_max_conta
            )
        ]
    )

    # Log information about high-quality bins
    thresholds = (
        f"(completeness >= {hq_min_completeness} and contamination <= {hq_max_conta})"
    )
    logger.info(
        f"{hq_bins}/{len(selected_bins)} selected bins have high quality {thresholds}"
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
            help="Resume mode: reuse existing temporary files if possible.",
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

    # Validate that exactly one of bin_dirs or contig2bin_tables is provided
    if bin_dirs is None and contig2bin_tables is None:
        raise typer.BadParameter(
            "Error: Either --bin-dirs or --contig2bin_tables must be provided. None were given."
        )

    if bin_dirs is not None and contig2bin_tables is not None:
        raise typer.BadParameter(
            "Error: Either --bin-dirs or --contig2bin_tables must be provided, but not both."
        )

    # High quality threshold used just to log number of high quality bins.
    hq_max_conta = 5
    hq_min_completeness = 90

    # Temporary files #
    out_tmp_dir: Path = outdir / "temporary_files"
    os.makedirs(out_tmp_dir, exist_ok=True)

    use_existing_protein_file = False

    faa_file = out_tmp_dir / "assembly_proteins.faa.gz"

    diamond_result_file = out_tmp_dir / "diamond_result.tsv.gz"

    # Output files #
    final_bin_report: Path = outdir / "final_bins_quality_reports.tsv"
    original_bin_report_dir: Path = outdir / "input_bins_quality_reports"

    if resume:
        io.check_resume_file(faa_file, diamond_result_file)
        use_existing_protein_file = True

    (
        contig_key_to_original_bin,
        contigs_in_bins,
        contig_to_length,
        contig_to_index,
    ) = parse_input_files(
        bin_dirs,
        contig2bin_tables,
        contigs,
        fasta_extensions=set(fasta_extensions),
    )

    if debug:
        index_to_contig_file = outdir / "index_to_contig.tsv"
        logger.info(f"Writing index to contig mapping to '{index_to_contig_file}'")
        with open(index_to_contig_file, "w") as flout:
            flout.write("\n".join((f"{i}\t{c}" for i, c in enumerate(contigs_in_bins))))

    if proteins and not resume:
        logger.info(f"Using the provided protein sequences file '{proteins}'")
        use_existing_protein_file = True

        cds.filter_faa_file(
            contigs_in_bins,
            input_faa_file=proteins,
            filtered_faa_file=faa_file,
        )

    contig_name_to_kegg_counter, contig_name_to_genes, contig_to_coding_length = (
        manage_protein_alignement(
            faa_file=faa_file,
            contigs_fasta=contigs,
            contigs_in_bins=contigs_in_bins,
            diamond_result_file=diamond_result_file,
            checkm2_db=checkm2_db,
            threads=threads,
            use_existing_protein_file=use_existing_protein_file,
            resume_diamond=resume,
            low_mem=low_mem,
        )
    )

    contig_to_kegg_counter = contig_manager.apply_contig_index(
        contig_to_index, contig_name_to_kegg_counter
    )
    contig_to_genes = contig_manager.apply_contig_index(
        contig_to_index, contig_name_to_genes
    )
    if contig_to_coding_length:
        contig_to_coding_length = contig_manager.apply_contig_index(
            contig_to_index, contig_to_coding_length
        )
    # Extract cds metadata ##
    logger.info("Computing CDS metadata")
    contig_metadat = cds.get_contig_cds_metadata(contig_to_genes, threads)

    contig_metadat["contig_to_kegg_counter"] = contig_to_kegg_counter
    contig_metadat["contig_to_length"] = contig_to_length

    logger.info("Adding size and assessing quality of input bins")
    original_bins = bin_quality.add_bin_metrics(
        list(contig_key_to_original_bin.values()),
        contig_metadat,
        contamination_weight,
        threads,
        disable_progress_bar=not progress or quiet,
    )
    contig_key_to_original_bin = {b.contigs_key: b for b in original_bins}

    bin_quality.add_bin_size_and_N50(original_bins, contig_to_length)

    if contig_to_coding_length:
        bin_quality.add_bin_coding_density(original_bins, contig_to_coding_length)

    logger.info(
        f"Writing original input bin metrics to directory '{original_bin_report_dir}'"
    )
    io.write_original_bin_metrics(original_bins, original_bin_report_dir)

    logger.info("Creating intermediate bins")

    contig_lengths = bin_quality.prepare_contig_sizes(contig_to_length)

    contig_key_to_new_bin = bin_manager.create_intermediate_bins(
        contig_key_to_original_bin,
        contig_lengths=contig_lengths,
        min_comp=min_completeness,
        max_conta=max_contamination,
        min_len=min_length,
        max_len=max_length,
        disable_progress_bar=not progress or quiet,
    )

    logger.info(f"Assessing quality for {len(contig_key_to_new_bin)} intermediate bins")

    new_bins = bin_quality.add_bin_metrics(
        bins=contig_key_to_new_bin.values(),
        contig_info=contig_metadat,
        contamination_weight=contamination_weight,
        threads=threads,
        disable_progress_bar=not progress or quiet,
    )
    contig_key_to_new_bin = {b.contigs_key: b for b in new_bins}

    contig_key_to_all_bin = contig_key_to_original_bin | contig_key_to_new_bin

    bin_quality.add_bin_size_and_N50(contig_key_to_all_bin.values(), contig_to_length)

    if debug:
        all_bin_compo_file = outdir / "all_bins_quality_reports.tsv"
        logger.info(f"Writing all bins to '{all_bin_compo_file}'")
        io.write_bin_info(
            contig_key_to_all_bin.values(), all_bin_compo_file, add_contigs=True
        )

    selected_bins = bin_manager.select_best_bins(
        contig_key_to_all_bin,
        min_completeness=min_completeness,
        max_contamination=max_contamination,
        prefix=prefix,
    )

    if contig_to_coding_length:
        bin_quality.add_bin_coding_density(selected_bins, contig_to_coding_length)

    logger.info(f"Writing selected bins information to '{final_bin_report}'")
    io.write_bin_info(selected_bins, output=final_bin_report)

    io.write_contig2bin_table(
        selected_bins,
        outdir / "final_contig_to_bin.tsv",
        contigs_in_bins,
    )

    if write_fasta_bins:
        logger.info(f"Writing selected bins FASTA files to '{outdir / 'final_bins'}'")
        io.write_bins_fasta(
            selected_bins,
            contigs,
            outdir=outdir / "final_bins",
            contigs_names=contigs_in_bins,
        )

    log_selected_bin_info(selected_bins, hq_min_completeness, hq_max_conta)

    return 0


def main():
    preprocess_args()

    app()
