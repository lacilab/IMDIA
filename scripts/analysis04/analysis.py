#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import deque
from dataclasses import replace
from pathlib import Path
import re
import subprocess
import sys
import yaml
from scripts.utils import composition
from scripts.utils.helpers import setup_logger

# A line written by a pipeline logger ("<date> <time>,<ms> - LEVEL - message").
LOG_LINE = re.compile(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3} - "
                      r"(DEBUG|INFO|WARNING|ERROR|CRITICAL) - ")
# Number of Snakemake output lines repeated in the error message when it fails.
TAIL_LINES = 30

def main(logfile: Path, datadir: Path, configfile: Path, pdbpath: Path, postpdb: Path,
         cores: int = None, rules: list[str] = None, unlock: bool = False,
         flags: dict | None = None):
    """
    Analyse MD results using a Snakemake-based workflow.

    Arguments:
        logfile    : Path to the log file.
        datadir    : Path to the directory containing simulation data.
        pdbpath    : Path to the starting PDB structure.
        postpdb    : Path to a PDB written from the simulation (chain IDs and numbering
                     as used in the analysis).
        configfile : Path to the YAML parameter file; its settings are copied into the
                     composition record, which is the Snakemake config file.
        cores      : Number of CPU cores to use.
        rules      : List of Snakemake rules to execute.
        unlock     : Unlock Snakemake workflow if it is locked.
        flags      : Resolved molecule-type flags (dna/rna/phos), stored in place of
                     the parameter file's values.
    """

    # Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug(f"Starting analysis with snakemake using config file: {configfile}.")

    ##### Map new chain IDs and numbering to the input structure ######
    # The composition record of the input PDB is kept next to it.
    record = Path(pdbpath).parent / composition.RECORD_NAME
    _record_chain_mapping(record, Path(postpdb), logger)
    # The analysis settings go into the same record, which is Snakemake's config file.
    _record_params(record, Path(configfile), flags, logger)

    if rules is not None:
        logger.debug(f"Executing specific rules: {rules}")
    else:
        logger.debug("No rules were specified. Executing all rules in the Snakefile.")

    # Run-specific values, passed on the command line next to the config file.
    config_args = [f"logfile={str(logfile)}", f"datadir={str(datadir)}"]

    if rules:
        config_args.append(f"rule_groups={' '.join(rules)}")

    if unlock:
        logger.debug("Unlocking snakemake workflow.")
        cmd = [
            "snakemake",
            "--snakefile", "scripts/analysis04/Snakefile",
            "--cores", "1",
            "--configfile", str(record),
            "--config",
            *config_args,
            "--unlock"
        ]
    else:
        cmd = [
            "snakemake",
            "--snakefile", "scripts/analysis04/Snakefile",
            "--cores", str(cores) if cores is not None else "all",
            "--configfile", str(record),
            "--config",
            *config_args
        ]

    _run_snakemake(cmd, logger)


def _run_snakemake(cmd: list[str], logger):
    """Run Snakemake and route its output.

    Lines written by the loggers of the rule scripts are printed to the console as they
    are: those loggers already chose which levels to show and wrote them to their own
    log files. Every other line (Snakemake's own output) is logged at DEBUG, so it is
    always in the log file and on the console only in verbose mode. Continuation lines
    of multi-line log messages and tracebacks count as Snakemake output. If Snakemake
    fails, its last lines are logged as an error.
    """
    # stderr is merged into stdout so the lines keep their order.
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    tail = deque(maxlen=TAIL_LINES)
    for line in process.stdout:
        line = line.rstrip("\n")
        tail.append(line)
        if LOG_LINE.match(line):
            # Rule-script log line: shown on the console only.
            print(line, flush=True)
        else:
            logger.debug(f"snakemake: {line}")

    # A non-zero exit code stops the stage, with the last lines as the reason.
    if process.wait() != 0:
        logger.error(f"Snakemake failed (exit code {process.returncode}). Last output:\n"
                     + "\n".join(tail))
        sys.exit(1)
    logger.info("Snakemake workflow completed.")


def _record_chain_mapping(record: Path, postpdb: Path, logger) -> dict:
    """Add the input -> simulation chain mapping to the composition record.

    Chains of the post-simulation PDB are matched by order with the chains of the
    input structure and verified by sequence; residues without a chain ID are split
    by the input chain lengths. The mapping is stored in the record as
    ``chain_mapping = {"post_pdb", "post_sha256", "all_verified", "chains"}``, where
    "chains" is keyed by position (see ``composition.map_to_structure``). Chains that
    cannot be verified are kept with ``verified: false`` and reported as warnings.
    Returns the "chains" part.
    """
    if not record.is_file():
        logger.error(f"Composition file {record} not found.")
        sys.exit(1)
    old = composition.load(record)

    # Pair the chains by order; a structure that cannot be read stops the stage.
    try:
        chains = composition.map_to_structure(old, postpdb)
    except Exception as err:
        logger.error(f"Could not read the post-simulation PDB {postpdb}: {err}")
        sys.exit(1)

    # Residues without a chain ID were assigned to chains by the input chain lengths.
    unnamed = [pos for pos, e in chains.items() if e["id_new"] is not None and not e["id_new"].strip()]
    if unnamed:
        logger.debug(f"{postpdb.name}: residues without a chain ID were split by the input "
                    f"chain lengths (0-based chain positions {', '.join(unnamed)}); select these "
                    "chains by their residue or atom ranges.")

    # One message per chain: confirmation, or a warning with the reason.
    for pos, e in chains.items():
        if e["verified"]:
            logger.debug(f"Chain {e['id_old']} (input) -> '{e['id_new']}' (simulation, "
                         f"0-based position {pos}, residues {e['resids_new'][0]}-{e['resids_new'][1]}), "
                         f"offset {e['offset']}")
        else:
            logger.warning(f"Chain at 0-based position {pos} ('{e['id_new']}' in {postpdb.name}) "
                           f"could not be mapped to the input structure: {e['note']}. "
                           "Simulation labels are used for this chain.")

    # Store the mapping with the fingerprint of the file it was derived from.
    mapping = {"post_pdb": str(postpdb), "post_sha256": composition.file_sha256(postpdb),
               "all_verified": all(e["verified"] for e in chains.values()),
               "chains": chains}
    composition.save(replace(old, chain_mapping=mapping), record)
    logger.debug(f"Chain mapping written to {record}")
    return chains


def _record_params(record: Path, configfile: Path, flags: dict | None, logger) -> dict:
    """Copy the parameter file into the composition record as ``params``.

    ``flags`` (protein/dna/rna/phos as resolved by the pipeline) replace the values of the
    parameter file, which may be "auto" or missing. ``structure`` and ``muts`` default
    to the folders of the input PDB (data/<structure>/<mut>/), and a single name is
    turned into a list. Returns the stored settings.
    """
    with open(configfile) as f:
        params = yaml.safe_load(f) or {}
    # Resolved molecule-type flags take precedence over the file's values.
    for key, value in (flags or {}).items():
        if params.get(key) != value:
            logger.debug(f"{key}: {value} (resolved; {Path(configfile).name} has "
                         f"{params.get(key)!r})")
        params[key] = value

    # The record JSON is in the mutation folder. Its parent is the structure folder.
    mutdir = Path(record).parent
    for key, default in (("structure", mutdir.parent.name), ("muts", mutdir.name)):
        value = params.get(key)
        # Missing or empty in params file: use the folder name.
        if not value:
            params[key] = [default]
            logger.debug(f"{key}: {params[key]} (from the PDB path)")
        elif isinstance(value, str):
            # A single name is stored as a one-element list.
            params[key] = [value]

    # Store the settings next to the composition and the chain mapping.
    comp = composition.load(record)
    composition.save(replace(comp, params=params), record)
    logger.debug(f"Analysis settings from {configfile} written to {record}.")
    return params


if __name__ == "__main__":
    main()
