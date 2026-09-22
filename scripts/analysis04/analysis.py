#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import shutil
import subprocess
import string
import MDAnalysis as mda
import yaml
from scripts.utils.helpers import setup_logger

def main(logfile: Path, datadir: Path, configfile: Path, pdbpath: Path, cores: int = None, rules: list[str] = None, unlock: bool = False):
    """
    Analyse MD results using a Snakemake-based workflow.

    Arguments:
        logfile    : Path to the log file.
        datadir    : Path to the directory containing simulation data.
        pdbpath    : Path to the starting PDB structure.
        configfile : Path to the YAML configuration file used by the Snakemake workflow.
        cores      : Number of CPU cores to use.
        rules      : List of Snakemake rules to execute.
        unlock     : Unlock Snakemake workflow if it is locked.
    """

    # Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug(f"Starting analysis with snakemake using config file: {configfile}.")

    ##### Get chain IDs and offsets ######
    # Extract chain ids from the PDB file
    uppercase_alphabet = list(string.ascii_uppercase)
    u = mda.Universe(pdbpath)
    chains = {}  # old_id -> {"id_new": ..., "offset": ...}
    i = 0

    for atom in u.atoms:
        if atom.segid:
            offset = int(atom.resid) - 1
            if atom.segid not in chains:
                new_id = uppercase_alphabet[i]
                chains[atom.segid] = {"id_new": new_id, "offset": offset}
                i += 1
                logger.debug(f"Found chain ID: '{atom.segid}' with starting residue number {atom.resid}. Assigned chain ID '{new_id}'.")
            elif offset < chains[atom.segid]['offset']:
                chains[atom.segid]['offset'] = offset
                logger.debug(f"Updated chain {atom.segid} (new ID {chains[atom.segid]['id_new']}) with new offset {offset}")

    # Final output keyed by new (alphabetical) ID
    chains_newid = {v["id_new"]: {"id_old": k, "offset": v["offset"]} for k, v in chains.items()}

    params_copy = shutil.copyfile(configfile, pdbpath.parent / "params_copy.yaml")

    with open(params_copy, "r+") as f:
        content = f.read()

        # Trim trailing whitespace/newlines
        cleaned = content.rstrip()

        # Rewrite cleaned content
        f.seek(0)
        f.write(cleaned)
        f.truncate()

        # Move to end and append new content
        f.write("\n\n")
        yaml.safe_dump({"chains": chains_newid}, f)
    
    if rules is not None:
        logger.debug(f"Executing specific rules: {rules}")
    else:
        logger.debug("No rules were specified. Executing all rules in the Snakefile.")
        
    config_args = [f"logfile={str(logfile)}", f"datadir={str(datadir)}"]

    if rules:
        config_args.append(f"rule_groups={' '.join(rules)}")

    if unlock:
        logger.debug("Unlocking snakemake workflow.")
        cmd = [
            "snakemake",
            "--snakefile", "scripts/analysis04/Snakefile",
            "--cores", "1",
            "--configfile", str(params_copy),    
            "--config",
            *config_args,
            "--unlock"
        ]
    else:
        cmd = [
            "snakemake",
            "--snakefile", "scripts/analysis04/Snakefile",
            "--cores", str(cores) if cores is not None else "all",
            "--configfile", str(params_copy),
            "--config",
            *config_args
        ]

    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()