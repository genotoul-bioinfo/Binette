"""
Module      : checkpoint
Description : Track intermediate results so that a run can be safely resumed.

Binette stores a small JSON manifest in its temporary directory. The manifest
records the identity of the run (fingerprints of its data inputs) and, for each
completed step, the inputs it was computed from and the files it produced
(size + sha256). A step is only reused when its inputs match the current run
and its output files are intact on disk.

Copyright   : (c) Jean Mainguy
License     : GPL-3.0
"""

import hashlib
import json
import logging
import os
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import binette as binette_init

logger = logging.getLogger(__name__)

FORMAT_VERSION = 1
CHECKPOINT_FILENAME = "checkpoint.json"

# Files and directories written by Binette in its output directory. Their presence
# without a checkpoint manifest indicates a run from an older Binette version.
BINETTE_OUTPUT_MARKERS = (
    "final_bins_quality_reports.tsv",
    "final_contig_to_bin.tsv",
    "final_bins",
    "input_bins_quality_reports",
    "temporary_files/assembly_proteins.faa.gz",
    "temporary_files/diamond_result.tsv.gz",
)

RUN_INPUT_LABELS = {
    "contigs": "contig sequences (--contigs)",
    "bins": "input bin sets (--bin_dirs / --contig2bin_tables)",
    "proteins": "protein file (--proteins)",
    "checkm2_db": "CheckM2 database (--checkm2_db)",
}


class OutputDirConflictError(Exception):
    """Raised when the output directory holds results from a different run."""


def sha256_file(path: Path, chunk_size: int = 1 << 20) -> str:
    """Compute the sha256 hex digest of a file."""
    digest = hashlib.sha256()
    with open(path, "rb") as fl:
        while chunk := fl.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def hash_json(obj: Any) -> str:
    """Compute a stable sha256 hex digest of a JSON-serialisable object."""
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def contigs_fingerprint(contig_to_seq_digest: dict[str, str]) -> str:
    """
    Fingerprint a set of contigs from per-contig sequence digests.

    The fingerprint is independent of the order of contigs in the FASTA file.

    :param contig_to_seq_digest: Mapping of contig name to the digest of its sequence.
    """
    digest = hashlib.sha256()
    for name in sorted(contig_to_seq_digest):
        digest.update(f"{name}\t{contig_to_seq_digest[name]}\n".encode())
    return digest.hexdigest()


def bins_fingerprint(bin_set_name_to_bins_info: dict[str, list[dict]]) -> str:
    """
    Fingerprint input bin sets: bin set names, bin names and their contigs.

    :param bin_set_name_to_bins_info: Mapping of bin set name to a list of bin info dicts
        with keys "bin_name" and "contigs".
    """
    normalized = {
        set_name: sorted(
            (str(info["bin_name"]), sorted(info["contigs"])) for info in bins_info
        )
        for set_name, bins_info in bin_set_name_to_bins_info.items()
    }
    return hash_json(normalized)


def checkm2_db_identity(db_path: Path) -> dict[str, Any]:
    """
    Identify a CheckM2 diamond database by file name and size.

    The database is too large to hash. The path is not used, so a database copied
    to another location (e.g. node-local scratch) is still recognised.
    """
    return {"name": db_path.name, "size": db_path.stat().st_size}


def partial_path(path: Path) -> Path:
    """Temporary path used to write a file before atomically moving it into place."""
    return path.with_name(f"partial_{path.name}")


def atomic_replace(tmp_path: Path, final_path: Path) -> None:
    """Move a fully written temporary file to its final location."""
    os.replace(tmp_path, final_path)


class Checkpoint:
    """Manifest of the intermediate results of a Binette run."""

    def __init__(self, manifest_path: Path, root: Path, data: dict[str, Any]):
        """
        :param manifest_path: Path of the JSON manifest.
        :param root: Directory that recorded output paths are relative to.
        :param data: Manifest content.
        """
        self.manifest_path = manifest_path
        self.root = root
        self.data = data

    @classmethod
    def new(
        cls, manifest_path: Path, root: Path, run_inputs: dict[str, Any]
    ) -> "Checkpoint":
        return cls(
            manifest_path,
            root,
            {
                "format_version": FORMAT_VERSION,
                "binette_version": binette_init.__version__,
                "run_inputs": run_inputs,
                "steps": {},
            },
        )

    @classmethod
    def load(cls, manifest_path: Path, root: Path) -> tuple["Checkpoint | None", str]:
        """
        Load a manifest from disk.

        :return: The checkpoint (or None) and a reason when it could not be loaded.
        """
        if not manifest_path.exists():
            return None, f"no checkpoint file '{manifest_path}'"
        try:
            with open(manifest_path) as fl:
                data = json.load(fl)
        except (OSError, json.JSONDecodeError) as exc:
            return None, f"checkpoint file '{manifest_path}' is unreadable ({exc})"

        if not isinstance(data, dict) or data.get("format_version") != FORMAT_VERSION:
            return None, (
                f"checkpoint file '{manifest_path}' has an unsupported format version"
            )
        data.setdefault("steps", {})
        return cls(manifest_path, root, data), ""

    @property
    def run_inputs(self) -> dict[str, Any]:
        return self.data.get("run_inputs", {})

    def save(self) -> None:
        """Write the manifest atomically."""
        self.manifest_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = partial_path(self.manifest_path)
        with open(tmp, "w") as fl:
            json.dump(self.data, fl, indent=2, sort_keys=True)
        atomic_replace(tmp, self.manifest_path)

    def _relative(self, path: Path) -> str:
        return Path(os.path.relpath(path, self.root)).as_posix()

    def step_status(self, step: str, inputs: dict[str, Any]) -> tuple[bool, str]:
        """
        Check whether a recorded step can be reused.

        :param step: Step name.
        :param inputs: Inputs of the step in the current run.
        :return: (is_valid, reason). The reason explains why a step is not reusable.
        """
        record = self.data["steps"].get(step)
        if record is None:
            return False, "no record of a completed previous run of this step"

        if record.get("inputs") != inputs:
            changed = sorted(
                key
                for key in set(inputs) | set(record.get("inputs", {}))
                if inputs.get(key) != record.get("inputs", {}).get(key)
            )
            return False, f"inputs changed since previous run ({', '.join(changed)})"

        for rel_path, expected in record.get("outputs", {}).items():
            path = self.root / rel_path
            if not path.is_file():
                return False, f"file '{path}' is missing"
            if path.stat().st_size != expected["size"]:
                return False, f"file '{path}' has changed or is incomplete"
            if sha256_file(path) != expected["sha256"]:
                return False, f"file '{path}' has changed"

        return True, ""

    def step_output_digest(self, step: str, path: Path) -> str | None:
        """Return the recorded sha256 of an output file of a step."""
        record = self.data["steps"].get(step, {})
        output = record.get("outputs", {}).get(self._relative(path))
        return None if output is None else output["sha256"]

    def record_step(
        self, step: str, inputs: dict[str, Any], outputs: Iterable[Path]
    ) -> None:
        """Record a completed step and save the manifest."""
        self.data["steps"][step] = {
            "inputs": inputs,
            "outputs": {
                self._relative(path): {
                    "size": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
                for path in outputs
            },
        }
        self.save()

    def invalidate(self, *steps: str) -> None:
        """Forget steps (e.g. before their outputs are rewritten) and save."""
        changed = False
        for step in steps:
            if self.data["steps"].pop(step, None) is not None:
                changed = True
        if changed:
            self.save()


def find_binette_outputs(outdir: Path) -> list[Path]:
    """List Binette output files/directories present in an output directory."""
    return [outdir / m for m in BINETTE_OUTPUT_MARKERS if (outdir / m).exists()]


def check_outdir(
    outdir: Path,
    manifest_path: Path,
    run_inputs: dict[str, Any],
    force: bool,
) -> Checkpoint:
    """
    Make sure the output directory can be used for the current run.

    - Empty/new output directory: a new checkpoint is created.
    - Output directory from a run with the same inputs: its checkpoint is returned.
    - Output directory from a run with different inputs, or from an older Binette
      version that cannot be verified: raises OutputDirConflictError, unless
      ``force`` is set, in which case a new checkpoint is created.

    :param outdir: Output directory of the run.
    :param manifest_path: Path of the checkpoint manifest.
    :param run_inputs: Fingerprints of the data inputs of the current run.
    :param force: Allow reusing an output directory from a different run.
    :raises OutputDirConflictError: When the output directory belongs to another run.
    """
    checkpoint, reason = Checkpoint.load(manifest_path, root=outdir)

    if checkpoint is not None:
        if checkpoint.run_inputs == run_inputs:
            logger.debug(f"Output directory '{outdir}' matches the current inputs")
            return checkpoint

        changed = [
            RUN_INPUT_LABELS.get(key, key)
            for key in sorted(set(run_inputs) | set(checkpoint.run_inputs))
            if run_inputs.get(key) != checkpoint.run_inputs.get(key)
        ]
        problem = (
            f"Output directory '{outdir}' contains results from a previous Binette run "
            f"with different inputs. Changed: {'; '.join(changed)}."
        )
    else:
        existing = find_binette_outputs(outdir)
        if not existing and not manifest_path.exists():
            checkpoint = Checkpoint.new(manifest_path, outdir, run_inputs)
            checkpoint.save()
            return checkpoint
        if manifest_path.exists():
            problem = (
                f"Output directory '{outdir}' contains results from a previous "
                f"Binette run that cannot be verified: {reason}."
            )
        else:
            problem = (
                f"Output directory '{outdir}' contains results from a previous Binette "
                "run that cannot be verified (no checkpoint information, "
                "probably produced by an older Binette version)."
            )

    if not force:
        raise OutputDirConflictError(
            f"{problem} Use a different output directory, or use --force to discard "
            "the previous intermediate files and overwrite the results."
        )

    logger.warning(
        f"{problem} --force is set: previous intermediate files are discarded "
        "and results will be overwritten."
    )
    checkpoint = Checkpoint.new(manifest_path, outdir, run_inputs)
    checkpoint.save()
    return checkpoint
