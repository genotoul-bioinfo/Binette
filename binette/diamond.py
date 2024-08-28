import subprocess
import logging
import sys
import shutil
import re
import pandas as pd
from collections import Counter

from checkm2 import keggData


def get_checkm2_db() -> str:
    """
    Get the path to the CheckM2 database.

    :return: The path to the CheckM2 database.
    """
    if shutil.which("checkm2") is None:
        logging.error("Make sure checkm2 is on your system path.")
        sys.exit(1)

    checkm2_database_raw = subprocess.run(["checkm2", "database", "--current"], text=True, stderr=subprocess.PIPE)

    if checkm2_database_raw.returncode != 0:
        logging.error(f"Something went wrong with checkm2:\n=======\n{checkm2_database_raw.stderr}========")
        sys.exit(1)

    reg_result = re.search("INFO: (/.*.dmnd)", checkm2_database_raw.stderr)

    if reg_result is None:
        logging.error(f"Something went wrong when retrieving checkm2 db path:\n{checkm2_database_raw.stderr}")
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
        raise FileNotFoundError(f"The '{tool_name}' tool is not found on your system PATH.")


def run(
    faa_file: str, output: str, db: str, log: str, threads: int = 1, query_cover: int = 80, subject_cover: int = 80,
    percent_id: int = 30, evalue: float = 1e-05, low_mem: bool = False):
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

    blocksize = 0.5 if low_mem else 2

    cmd = (
        "diamond blastp --outfmt 6 --max-target-seqs 1 "
        f"--query {faa_file} "
        f"-o {output} "
        f"--threads {threads} "
        f"--db {db} "
        f"--query-cover {query_cover} "
        f"--subject-cover {subject_cover} "
        f"--id {percent_id} "
        f"--evalue {evalue} --block-size {blocksize} 2> {log}"
    )

    logging.info("Running diamond")
    logging.info(cmd)

    run = subprocess.run(cmd, shell=True)

    if run.returncode != 0:
        logging.error(f"An error occurred while running DIAMOND. Check log file: {log}")
        sys.exit(1)

    logging.info("Finished Running DIAMOND")

def get_contig_to_kegg_id(diamond_result_file: str) -> dict:
    """
    Get a dictionary mapping contig IDs to KEGG annotations from a Diamond result file.

    :param diamond_result_file: Path to the Diamond result file.
    :return: A dictionary mapping contig IDs to KEGG annotations.
    """
    diamon_results_df = pd.read_csv(diamond_result_file, sep="\t", usecols=[0, 1], names=["ProteinID", "annotation"])
    diamon_results_df[["Ref100_hit", "Kegg_annotation"]] = diamon_results_df["annotation"].str.split(
        "~", n=1, expand=True
    )

    KeggCalc = keggData.KeggCalculator()
    defaultKOs = KeggCalc.return_default_values_from_category("KO_Genes")

    diamon_results_df = diamon_results_df.loc[diamon_results_df["Kegg_annotation"].isin(defaultKOs.keys())]
    diamon_results_df["contig"] = diamon_results_df["ProteinID"].str.split("_", n=-1).str[:-1].str.join("_")

    contig_to_kegg_counter = (
        diamon_results_df.groupby("contig").agg({"Kegg_annotation": Counter}).reset_index()
    )

    contig_to_kegg_counter = dict(zip(contig_to_kegg_counter["contig"], contig_to_kegg_counter["Kegg_annotation"]))

    return contig_to_kegg_counter
