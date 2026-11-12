import inspect
import logging
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from checkm2 import keggData

from binette.checkpoint import atomic_replace, partial_path

logger = logging.getLogger(__name__)


def get_checkm2_db() -> str:
    """
    Get the path to the CheckM2 database.

    :return: The path to the CheckM2 database.
    """
    if shutil.which("checkm2") is None:
        logger.error("Make sure checkm2 is on your system path")
        sys.exit(1)

    checkm2_database_raw = subprocess.run(
        ["checkm2", "database", "--current"], text=True, stderr=subprocess.PIPE
    )

    if checkm2_database_raw.returncode != 0:
        logger.error(
            f"Something went wrong with checkm2:\n=======\n{checkm2_database_raw.stderr}========"
        )
        sys.exit(1)

    reg_result = re.search("INFO: (/.*.dmnd)", checkm2_database_raw.stderr)

    if reg_result is None:
        logger.error(
            f"Something went wrong when retrieving checkm2 db path:\n{checkm2_database_raw.stderr}"
        )
        sys.exit(1)
    else:
        db_path = reg_result.group(1)

    return db_path


def check_tool_exists(tool_name: str):
    """
    Check if a specified tool is on the system's PATH.

    :param tool_name: The name of the tool to check for.
    :type tool_name: str
    :raises FileNotFoundError: If the tool is not found on the system's PATH.
    """
    if shutil.which(tool_name) is None:
        raise FileNotFoundError(
            f"The '{tool_name}' tool is not found on your system PATH."
        )


def resolve_checkm2_db(checkm2_db: Path | None) -> Path:
    """
    Return the path of the CheckM2 diamond database to use.

    :param checkm2_db: Database given by the user, or None to use the one set via `checkm2 database`.
    :raises FileNotFoundError: If the given database does not exist.
    """
    if checkm2_db is None:
        return Path(get_checkm2_db())
    if not checkm2_db.exists():
        raise FileNotFoundError(checkm2_db)
    return checkm2_db


def get_block_size(low_mem: bool) -> float:
    """
    DIAMOND block size used by Binette. It changes memory usage and can also
    slightly change alignment results.

    :param low_mem: Use low memory mode if True.
    """
    return 0.5 if low_mem else 2


def run(
    faa_file: str,
    output: str,
    db: str,
    log: str,
    threads: int = 1,
    query_cover: int = 80,
    subject_cover: int = 80,
    percent_id: int = 30,
    evalue: float = 1e-05,
    low_mem: bool = False,
):
    """
    Run Diamond with specified parameters.

    :param faa_file: Path to the input protein sequence file (FASTA format).
    :param output: Path to the Diamond output file.
    :param db: Path to the Diamond database.
    :param log: Path to the log file.
    :param threads: Number of CPU threads to use (default is 1).
    :param query_cover: Minimum query coverage percentage (default is 80).
    :param subject_cover: Minimum subject coverage percentage (default is 80).
    :param percent_id: Minimum percent identity (default is 30).
    :param evalue: Maximum e-value threshold (default is 1e-05).
    :param low_mem: Use low memory mode if True (default is False).
    """
    check_tool_exists("diamond")

    # DIAMOND writes to a temporary file which is moved to its final location
    # only on success, so an interrupted job never leaves a truncated result file.
    final_output = Path(output)
    tmp_output = partial_path(final_output)

    blocksize = get_block_size(low_mem)

    cmd = (
        "diamond blastp --outfmt 6 --max-target-seqs 1 "
        f"--query {faa_file} "
        f"-o {tmp_output} "
        f"--threads {threads} "
        f"--db {db} "
        f"--compress 1 "
        f"--query-cover {query_cover} "
        f"--subject-cover {subject_cover} "
        f"--id {percent_id} "
        f"--evalue {evalue} --block-size {blocksize} 2> {log}"
    )

    logger.info("Running diamond")
    logger.info(f"Command: {cmd}")

    run = subprocess.run(cmd, shell=True)

    if run.returncode != 0:
        logger.error(f"An error occurred while running DIAMOND. Check log file '{log}'")
        sys.exit(1)

    atomic_replace(tmp_output, final_output)
    logger.info("Finished running DIAMOND")


def get_contig_to_kegg_id(diamond_result_file: str) -> dict:
    """
    Get a dictionary mapping contig IDs to KEGG annotations from a Diamond result file.

    :param diamond_result_file: Path to the Diamond result file.
    :return: A dictionary mapping contig IDs to KEGG annotations.
    """
    diamond_results_df = pd.read_csv(
        diamond_result_file, sep="\t", usecols=[0, 1], names=["ProteinID", "annotation"]
    )

    if diamond_results_df.empty:
        logger.error(
            f"DIAMOND result file '{diamond_result_file}' is empty. "
            "This can happen with low-quality assemblies where DIAMOND produces no hits."
        )
        sys.exit(3)

    diamond_results_df[["Ref100_hit", "Kegg_annotation"]] = diamond_results_df[
        "annotation"
    ].str.split("~", n=1, expand=True)

    KeggCalc = keggData.KeggCalculator()
    defaultKOs = KeggCalc.return_default_values_from_category("KO_Genes")

    diamond_results_df = diamond_results_df.loc[
        diamond_results_df["Kegg_annotation"].isin(defaultKOs.keys())
    ]
    diamond_results_df["contig"] = (
        diamond_results_df["ProteinID"].str.split("_", n=-1).str[:-1].str.join("_")
    )

    contig_to_kegg_counter = (
        diamond_results_df.groupby("contig")
        .agg({"Kegg_annotation": Counter})
        .reset_index()
    )

    contig_to_kegg_counter = dict(
        zip(
            contig_to_kegg_counter["contig"],
            contig_to_kegg_counter["Kegg_annotation"],
            strict=False,
        )
    )

    return contig_to_kegg_counter


def result_parameters(low_mem: bool) -> dict:
    """
    DIAMOND parameters that can change its results. Threads do not, but the block
    size set by --low_mem does.

    :param low_mem: Whether DIAMOND low memory mode is used.
    """
    params = inspect.signature(run).parameters
    result_params = {
        name: params[name].default
        for name in ("query_cover", "subject_cover", "percent_id", "evalue")
    }
    result_params["block_size"] = get_block_size(low_mem)
    return result_params
