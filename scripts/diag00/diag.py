#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from pathlib import Path
from scripts.utils.helpers import setup_logger
from scripts.diag00 import diag_01_checkgaps as checkgaps


def main(logfile: Path, datadir: Path, filename: str, stop_on_gaps: bool = True):
    """
    Diagnose a PDB structure before preparing it for MD simulations.

    Steps:
        1.Gap check: chain breaks, backbone distortions, residues with missing heavy atoms
          and residue-numbering jumps (written to diag/<file>_gaps.csv).

    Arguments:
        logfile      : Path to the log file.
        datadir      : Directory containing the input PDB file.
        filename     : Name of the PDB file (with or without .pdb)
        stop_on_gaps : Stop the pipeline if chain breaks or backbone distortions are found
                       (default: True). Missing atoms never stop it.
    """

    ##### INITIALIZE LOG
    logger = setup_logger("diag", logfile)

    ###### CHECK INPUTS ######
    # Normalize filename: accept with or without .pdb, store basename without extension
    file = Path(filename).stem

    # check datadir exists
    if not datadir.is_dir():
        logger.error(f"datadir '{datadir}' does not exist or is not a directory.")
        sys.exit(1)

    # check input PDB file exists
    path_pdb = datadir / f"{file}.pdb"
    if not path_pdb.is_file():
        logger.error(f"Input PDB file '{str(path_pdb)}' not found.")
        sys.exit(1)

    logger.info(f"Diagnosing {datadir}/{file}.pdb")

    # Output folder for the diagnosis files.
    path_diag = datadir / "diag"
    if path_diag.is_dir():
        logger.debug(f"{str(path_diag)} already exists.")
    else:
        logger.debug(f"{str(path_diag)} does not exist. Creating new directory.")
        path_diag.mkdir(exist_ok=True)

    ##### DIAGNOSIS ######

    ### 1. Gap check
    checkgaps.main(logfile=logfile, pdb=path_pdb, outdir=path_diag, stop_on_gaps=stop_on_gaps)

    logger.info("Diagnosis completed.")


if __name__ == "__main__":
    main()
