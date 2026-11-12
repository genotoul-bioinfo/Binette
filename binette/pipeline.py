"""
Module      : pipeline
Description : Orchestration of a Binette run, from input parsing to final bin selection.
Copyright   : (c) Jean Mainguy
License     : GPL-3.0
"""

import hashlib
import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pyfastx

import binette as binette_init
from binette import bin_manager, bin_quality, cds, checkpoint, contig_manager, diamond
from binette import io_manager as io
from binette.bin_manager import Bin
from binette.config import OutputPaths, RunConfig

logger = logging.getLogger(__name__)

# High quality thresholds, only used to log the number of high quality bins.
HQ_MIN_COMPLETENESS = 90
HQ_MAX_CONTAMINATION = 5


@dataclass
class InputData:
    """Bins and contigs parsed from the input files."""

    contig_key_to_original_bin: dict[Any, Bin]
    contigs_in_bins: list[str]
    contig_to_length: dict[int, int]
    contig_to_index: dict[str, int]
    fingerprints: dict[str, str]


def run(config: RunConfig) -> int:
    """
    Run Binette.

    :param config: Parameters of the run.
    :raises checkpoint.OutputDirConflictError: If the output directory holds results
        of a run with different inputs (and ``config.force`` is not set).
    :return: Exit code.
    """
    paths = OutputPaths.from_outdir(config.outdir)

    inputs = parse_input_files(
        config.bin_dirs,
        config.contig2bin_tables,
        config.contigs,
        fasta_extensions=set(config.fasta_extensions),
    )

    diamond_db_path = diamond.resolve_checkm2_db(config.checkm2_db)
    run_inputs = make_run_inputs(config, inputs, diamond_db_path)

    ckpt = checkpoint.check_outdir(
        config.outdir, paths.manifest, run_inputs, config.force
    )

    if config.resume and not ckpt.data["steps"]:
        logger.info(
            "No intermediate results from a previous run to resume from: starting from scratch"
        )

    if is_already_completed(config, paths, ckpt, run_inputs):
        return 0

    # Results of this run will be rewritten: they no longer reflect a completed run.
    ckpt.invalidate("final", "bins_fasta")

    if config.debug:
        write_index_to_contig(inputs.contigs_in_bins, paths.index_to_contig_file)

    contig_metadata, contig_to_coding_length, diamond_step_inputs = annotate_contigs(
        config, paths, inputs, ckpt, run_inputs, diamond_db_path
    )

    contig_key_to_original_bin = assess_input_bins(
        config, paths, inputs, contig_metadata, contig_to_coding_length
    )

    contig_key_to_all_bin = create_and_assess_intermediate_bins(
        config, paths, inputs, contig_key_to_original_bin, contig_metadata
    )

    selected_bins = select_bins(config, contig_key_to_all_bin, contig_to_coding_length)

    write_final_outputs(
        config,
        paths,
        inputs,
        selected_bins,
        ckpt,
        final_step_inputs(config, diamond_step_inputs),
    )

    return 0


# Inputs #


def parse_input_files(
    bin_dirs: list[Path] | None,
    contig2bin_tables: list[Path] | None,
    contigs_fasta: Path,
    fasta_extensions: set[str] | None = None,
) -> InputData:
    """
    Parses input files to retrieve information related to bins and contigs.

    :param bin_dirs: List of paths to directories containing bin FASTA files.
    :param contig2bin_tables: List of paths to contig-to-bin tables.
    :param contigs_fasta: Path to the contigs FASTA file.
    :param fasta_extensions: Possible fasta extensions to look for in the bin directory.

    :return: Parsed bins and contigs, with input fingerprints ("contigs" and "bins")
        identifying the run.
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

    logger.info(
        f"Parsing contig fasta file '{contigs_fasta}' to retrieve contig lengths"
    )

    contigs_in_bins_set = set(contigs_in_bins)
    contig_to_length = {}
    contig_to_seq_digest = {}
    for name, seq in pyfastx.Fastx(contigs_fasta.as_posix()):
        if name in contigs_in_bins_set:
            contig_to_length[name] = len(seq)
            contig_to_seq_digest[name] = hashlib.sha256(seq.encode()).hexdigest()

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

    return InputData(
        contig_key_to_original_bin=contig_key_to_bin,
        contigs_in_bins=contigs_in_bins,
        contig_to_length=contig_id_to_length,
        contig_to_index=contig_to_index,
        fingerprints={
            "contigs": checkpoint.contigs_fingerprint(contig_to_seq_digest),
            "bins": checkpoint.bins_fingerprint(bin_set_name_to_bins_info),
        },
    )


def write_index_to_contig(contigs_in_bins: list[str], output_file: Path) -> None:
    """Write the mapping of contig index to contig name (debug mode)."""
    logger.info(f"Writing index to contig mapping to '{output_file}'")
    with open(output_file, "w") as flout:
        flout.write("\n".join((f"{i}\t{c}" for i, c in enumerate(contigs_in_bins))))


# Checkpoint steps #


def make_run_inputs(
    config: RunConfig, inputs: InputData, diamond_db_path: Path
) -> dict[str, Any]:
    """
    Identity of the run: data inputs only. Parameters that only affect bin
    selection can change between runs sharing an output directory.
    """
    return {
        "contigs": inputs.fingerprints["contigs"],
        "bins": inputs.fingerprints["bins"],
        "proteins": checkpoint.sha256_file(config.proteins)
        if config.proteins
        else None,
        "checkm2_db": checkpoint.checkm2_db_identity(diamond_db_path),
    }


def make_diamond_step_inputs(
    ckpt: checkpoint.Checkpoint, faa_file: Path, checkm2_db: dict, low_mem: bool
) -> dict:
    """
    Inputs of the DIAMOND step: the proteins it runs on (as recorded in the
    checkpoint), the CheckM2 database and the DIAMOND parameters.
    """
    return {
        "proteins_sha256": ckpt.step_output_digest("proteins", faa_file),
        "checkm2_db": checkm2_db,
        "parameters": diamond.result_parameters(low_mem),
    }


def final_step_inputs(config: RunConfig, diamond_step_inputs: dict) -> dict:
    """
    Inputs of the final results step. It chains to upstream steps: final results are
    only reusable if the DIAMOND results (and through them, the proteins) match.
    """
    return {
        "binette_version": binette_init.__version__,
        "parameters": config.selection_parameters(),
        "diamond": diamond_step_inputs,
    }


def bins_fasta_step_inputs(contig2bin_table: Path) -> dict:
    """Inputs of the bins FASTA step: the final bins it was written from."""
    return {"contig2bin_table_sha256": checkpoint.sha256_file(contig2bin_table)}


def record_bins_fasta_step(
    ckpt: checkpoint.Checkpoint, contig2bin_table: Path, final_bins_dir: Path
) -> None:
    """Record the final bins FASTA files in the checkpoint."""
    ckpt.record_step(
        "bins_fasta",
        bins_fasta_step_inputs(contig2bin_table),
        sorted(final_bins_dir.glob("*.fa")),
    )


def ensure_bins_fasta(
    ckpt: checkpoint.Checkpoint,
    contig2bin_table: Path,
    contigs_fasta: Path,
    final_bins_dir: Path,
) -> None:
    """
    Make sure the FASTA files of the final bins of a completed run are present and
    intact. Missing or modified files are written again from the final contig2bin
    table, without recomputing bin selection.

    :param ckpt: Checkpoint of the run.
    :param contig2bin_table: Final contig2bin table of the completed run.
    :param contigs_fasta: Path to the contigs FASTA file.
    :param final_bins_dir: Directory of the final bins FASTA files.
    """
    is_valid, reason = ckpt.step_status(
        "bins_fasta", bins_fasta_step_inputs(contig2bin_table)
    )
    if is_valid:
        return

    logger.info(
        f"Writing selected bins FASTA files to '{final_bins_dir}' from "
        f"'{contig2bin_table}' ({reason})"
    )
    for fasta_file in final_bins_dir.glob("*.fa"):
        fasta_file.unlink()
    io.write_bins_fasta_from_contig_to_bin(
        io.read_contig2bin_table(contig2bin_table), contigs_fasta, final_bins_dir
    )
    record_bins_fasta_step(ckpt, contig2bin_table, final_bins_dir)


def is_already_completed(
    config: RunConfig,
    paths: OutputPaths,
    ckpt: checkpoint.Checkpoint,
    run_inputs: dict[str, Any],
) -> bool:
    """
    Check whether a previous run completed with the same inputs and parameters.
    If so, missing bin FASTA files are written again when requested.
    """
    diamond_step_inputs = make_diamond_step_inputs(
        ckpt, paths.faa_file, run_inputs["checkm2_db"], config.low_mem
    )
    if not ckpt.can_reuse(
        "final",
        final_step_inputs(config, diamond_step_inputs),
        config.resume,
        "final results",
    ):
        return False

    if config.write_fasta_bins:
        ensure_bins_fasta(
            ckpt, paths.final_contig2bin_table, config.contigs, paths.final_bins_dir
        )
    logger.info(
        f"Binette already completed with identical inputs and parameters in '{config.outdir}': "
        "nothing to do. Use --no-resume to rerun."
    )
    return True


# Proteins and DIAMOND #


def get_proteins(
    faa_file: Path,
    coding_length_file: Path,
    contigs_fasta: Path,
    contigs_in_bins: set[str],
    proteins: Path | None,
    threads: int,
    reuse: bool,
) -> tuple[dict[str, list[str]], dict[str, int | None] | None, list[Path]]:
    """
    Predict genes, filter user-provided proteins, or reuse proteins from a previous run.

    :param faa_file: Path of the protein file used by Binette.
    :param coding_length_file: Path of the file storing contig coding lengths.
    :param contigs_fasta: Path to the contigs FASTA file.
    :param contigs_in_bins: Set of contigs in bins.
    :param proteins: Protein file provided by the user, or None to predict genes.
    :param threads: Number of threads for gene prediction.
    :param reuse: Reuse files of a previous run instead of computing them.

    :return: contig_to_genes, contig_to_coding_length (None when not computable) and
        the list of files produced.
    """
    if reuse:
        logger.info(f"Parsing protein sequences from '{faa_file}'")
        contig_to_genes = cds.parse_faa_file(faa_file.as_posix())
        contig_to_coding_len = (
            None if proteins else io.read_contig_coding_length(coding_length_file)
        )

    elif proteins:
        logger.info(f"Using the provided protein sequences file '{proteins}'")
        cds.filter_faa_file(
            contigs_in_bins,
            input_faa_file=proteins,
            filtered_faa_file=faa_file,
        )
        logger.info(f"Parsing protein sequences from '{faa_file}'")
        contig_to_genes = cds.parse_faa_file(faa_file.as_posix())
        contig_to_coding_len = None

    else:
        logger.info(f"Loading {len(contigs_in_bins)} contigs from '{contigs_fasta}'")
        contigs_iterator = (
            (name, seq)
            for name, seq in pyfastx.Fastx(contigs_fasta.as_posix())
            if name in contigs_in_bins
        )
        contig_to_genes, contig_to_coding_len = cds.predict(
            contigs_iterator, faa_file.as_posix(), threads
        )
        io.write_contig_coding_length(contig_to_coding_len, coding_length_file)

    if proteins:
        io.check_contig_consistency(
            contigs_in_bins,
            contig_to_genes,
            contigs_fasta.as_posix(),
            faa_file.as_posix(),
        )
        logger.info(
            "Coding density will not be computed (using provided protein sequences)"
        )
        outputs = [faa_file]
    else:
        logger.info("Coding density will be computed from identified genes")
        outputs = [faa_file, coding_length_file]

    return contig_to_genes, contig_to_coding_len, outputs


def get_diamond_hits(
    faa_file: Path,
    contigs_fasta: Path,
    contigs_in_bins: set[str],
    diamond_result_file: Path,
    diamond_db_path: Path,
    threads: int,
    low_mem: bool,
    reuse: bool,
) -> dict[str, Counter]:
    """
    Run DIAMOND on proteins against the CheckM2 database (or reuse a previous result)
    and parse its results.

    :param faa_file: Path to the protein file.
    :param contigs_fasta: Path to the contigs FASTA file.
    :param contigs_in_bins: Set of contigs in bins.
    :param diamond_result_file: Path of the DIAMOND result file.
    :param diamond_db_path: Path of the CheckM2 diamond database.
    :param threads: Number of threads.
    :param low_mem: Use DIAMOND low memory mode.
    :param reuse: Reuse the result file of a previous run instead of running DIAMOND.

    :return: Mapping of contig name to a counter of KEGG annotations.
    """
    if not reuse:
        diamond_log = (
            diamond_result_file.parents[0]
            / f"{diamond_result_file.stem.split('.')[0]}.log"
        )

        diamond.run(
            faa_file.as_posix(),
            diamond_result_file.as_posix(),
            diamond_db_path.as_posix(),
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

    return contig_to_kegg_counter


def annotate_contigs(
    config: RunConfig,
    paths: OutputPaths,
    inputs: InputData,
    ckpt: checkpoint.Checkpoint,
    run_inputs: dict[str, Any],
    diamond_db_path: Path,
) -> tuple[dict, dict[int, int | None] | None, dict]:
    """
    Get proteins and DIAMOND hits (computed or reused from a previous run) and build
    the contig metadata used to assess bin quality.

    :return: contig metadata, contig coding lengths (None when using user-provided
        proteins), both indexed by contig index, and the inputs of the DIAMOND step.
    """
    contigs_in_bins = set(inputs.contigs_in_bins)

    # Proteins #
    proteins_step_inputs = {
        "contigs": run_inputs["contigs"],
        "proteins": run_inputs["proteins"],
    }
    reuse_proteins = ckpt.can_reuse(
        "proteins", proteins_step_inputs, config.resume, "protein sequences"
    )
    if not reuse_proteins:
        ckpt.invalidate("proteins", "diamond")

    contig_name_to_genes, contig_to_coding_length, protein_outputs = get_proteins(
        faa_file=paths.faa_file,
        coding_length_file=paths.coding_length_file,
        contigs_fasta=config.contigs,
        contigs_in_bins=contigs_in_bins,
        proteins=config.proteins,
        threads=config.threads,
        reuse=reuse_proteins,
    )
    if not reuse_proteins:
        ckpt.record_step("proteins", proteins_step_inputs, protein_outputs)

    # DIAMOND #
    diamond_step_inputs = make_diamond_step_inputs(
        ckpt, paths.faa_file, run_inputs["checkm2_db"], config.low_mem
    )
    reuse_diamond = ckpt.can_reuse(
        "diamond", diamond_step_inputs, config.resume, "DIAMOND results"
    )
    if not reuse_diamond:
        ckpt.invalidate("diamond")

    contig_name_to_kegg_counter = get_diamond_hits(
        faa_file=paths.faa_file,
        contigs_fasta=config.contigs,
        contigs_in_bins=contigs_in_bins,
        diamond_result_file=paths.diamond_result_file,
        diamond_db_path=diamond_db_path,
        threads=config.threads,
        low_mem=config.low_mem,
        reuse=reuse_diamond,
    )
    if not reuse_diamond:
        ckpt.record_step("diamond", diamond_step_inputs, [paths.diamond_result_file])

    # Contig metadata #
    contig_to_index = inputs.contig_to_index
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

    logger.info("Computing CDS metadata")
    contig_metadata = cds.get_contig_cds_metadata(contig_to_genes, config.threads)
    contig_metadata["contig_to_kegg_counter"] = contig_to_kegg_counter
    contig_metadata["contig_to_length"] = inputs.contig_to_length

    return contig_metadata, contig_to_coding_length, diamond_step_inputs


# Bins #


def assess_input_bins(
    config: RunConfig,
    paths: OutputPaths,
    inputs: InputData,
    contig_metadata: dict,
    contig_to_coding_length: dict[int, int | None] | None,
) -> dict[Any, Bin]:
    """
    Assess quality of the input bins and write their metrics.

    :return: Mapping of contig key to input bins with quality metrics.
    """
    logger.info("Adding size and assessing quality of input bins")
    original_bins = bin_quality.add_bin_metrics(
        list(inputs.contig_key_to_original_bin.values()),
        contig_metadata,
        config.contamination_weight,
        config.threads,
        disable_progress_bar=not config.show_progress,
    )

    bin_quality.add_bin_size_and_N50(original_bins, inputs.contig_to_length)

    if contig_to_coding_length:
        bin_quality.add_bin_coding_density(original_bins, contig_to_coding_length)

    io.clear_previous_outputs(paths.final_bins_dir, paths.original_bin_report_dir)

    logger.info(
        f"Writing original input bin metrics to directory '{paths.original_bin_report_dir}'"
    )
    io.write_original_bin_metrics(original_bins, paths.original_bin_report_dir)

    return {b.contigs_key: b for b in original_bins}


def create_and_assess_intermediate_bins(
    config: RunConfig,
    paths: OutputPaths,
    inputs: InputData,
    contig_key_to_original_bin: dict[Any, Bin],
    contig_metadata: dict,
) -> dict[Any, Bin]:
    """
    Create intermediate bins from the input bins and assess their quality.

    :return: Mapping of contig key to all bins (input and intermediate).
    """
    logger.info("Creating intermediate bins")

    contig_lengths = bin_quality.prepare_contig_sizes(inputs.contig_to_length)

    contig_key_to_new_bin = bin_manager.create_intermediate_bins(
        contig_key_to_original_bin,
        contig_lengths=contig_lengths,
        min_comp=config.min_completeness,
        max_conta=config.max_contamination,
        min_len=config.min_length,
        max_len=config.max_length,
        disable_progress_bar=not config.show_progress,
    )

    logger.info(f"Assessing quality for {len(contig_key_to_new_bin)} intermediate bins")

    new_bins = bin_quality.add_bin_metrics(
        bins=contig_key_to_new_bin.values(),
        contig_info=contig_metadata,
        contamination_weight=config.contamination_weight,
        threads=config.threads,
        disable_progress_bar=not config.show_progress,
    )
    contig_key_to_new_bin = {b.contigs_key: b for b in new_bins}

    contig_key_to_all_bin = contig_key_to_original_bin | contig_key_to_new_bin

    bin_quality.add_bin_size_and_N50(
        contig_key_to_all_bin.values(), inputs.contig_to_length
    )

    if config.debug:
        logger.info(f"Writing all bins to '{paths.all_bin_compo_file}'")
        io.write_bin_info(
            contig_key_to_all_bin.values(), paths.all_bin_compo_file, add_contigs=True
        )

    return contig_key_to_all_bin


def select_bins(
    config: RunConfig,
    contig_key_to_all_bin: dict[Any, Bin],
    contig_to_coding_length: dict[int, int | None] | None,
) -> list[Bin]:
    """Select the best non-overlapping bins."""
    selected_bins = bin_manager.select_best_bins(
        contig_key_to_all_bin,
        min_completeness=config.min_completeness,
        max_contamination=config.max_contamination,
        prefix=config.prefix,
    )

    if contig_to_coding_length:
        bin_quality.add_bin_coding_density(selected_bins, contig_to_coding_length)

    return selected_bins


def write_final_outputs(
    config: RunConfig,
    paths: OutputPaths,
    inputs: InputData,
    selected_bins: list[Bin],
    ckpt: checkpoint.Checkpoint,
    final_inputs: dict,
) -> None:
    """Write the final results and record them in the checkpoint."""
    logger.info(f"Writing selected bins information to '{paths.final_bin_report}'")
    io.write_bin_info(selected_bins, output=paths.final_bin_report)

    io.write_contig2bin_table(
        selected_bins,
        paths.final_contig2bin_table,
        inputs.contigs_in_bins,
    )

    if config.write_fasta_bins:
        logger.info(f"Writing selected bins FASTA files to '{paths.final_bins_dir}'")
        io.write_bins_fasta(
            selected_bins,
            config.contigs,
            outdir=paths.final_bins_dir,
            contigs_names=inputs.contigs_in_bins,
        )

    bin_quality.log_selected_bin_info(
        selected_bins, HQ_MIN_COMPLETENESS, HQ_MAX_CONTAMINATION
    )

    final_outputs = [
        paths.final_bin_report,
        paths.final_contig2bin_table,
        *sorted(paths.original_bin_report_dir.glob("input_bins_*.tsv")),
    ]
    if config.debug:
        final_outputs += [paths.index_to_contig_file, paths.all_bin_compo_file]
    ckpt.record_step("final", final_inputs, final_outputs)

    if config.write_fasta_bins:
        record_bins_fasta_step(ckpt, paths.final_contig2bin_table, paths.final_bins_dir)
