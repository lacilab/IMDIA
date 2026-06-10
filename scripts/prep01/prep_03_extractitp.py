#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import re
import glob
import shutil
import time
from scripts.utils.helpers import setup_logger
from pathlib import Path


def main(logfile: Path, top: Path, dna: bool = False):

    """
    Generate separate .itp files from a GROMACS topology file.

    Arguments:
        logfile     : Path to the log file.
        top         : Path to the GROMACS topology file
        dna         : True if the system contains DNA.
    """

    gmxdir = top.parent
    path_toppar = gmxdir / "toppar"
    expected_itps = 8 if dna else 6

    ##### Initialize logger
    logger = setup_logger("prep", logfile)

    # Check for checkpoint files
    checkpoint_files = list(path_toppar.glob(f"check_03_extractitp_*.done"))
    if len(checkpoint_files) > 0:
        logger.warning(f"Checkpoint file(s) found: {', '.join([str(f) for f in checkpoint_files])}")
        itps = list(path_toppar.glob(f"*.itp"))
        if len(itps) != expected_itps:
            logger.warning(f"Expected {expected_itps} .itp files, but found {len(itps)}.")
        else:
            logger.info(f"The number of .itp files found matches the expected count of ({expected_itps}).")
        logger.warning("Skipping .itp extraction.")
        return

    # Check that gmxdir/toppar does not exist
    if not path_toppar.is_dir():
        path_toppar.mkdir(exist_ok=True)
        if path_toppar.is_dir():
            logger.info(f"{str(path_toppar)} directory created.")
        else:
            logger.error(f"{str(path_toppar)} directory could not be created.")
            sys.exit(1)
    else:
        logger.warning(f"{str(path_toppar)} already exists")

    # Extract molecule itps
    pattern = re.compile(r"\[\s*moleculetype\s*\]")
    with open(top, "r") as f:
        n = sum(1 for line in f if pattern.search(line))
    logger.debug(f"Found {n} molecule types in {top}.")
    i = 0
    temp_file = top.parent / f"{top.name}.tmp"
    inside = False
    dict_mols = {}
    with open(top, 'r') as fin, open(temp_file, "w") as fout:
        for line in fin:
            if "[ moleculetype ]" in line:
                inside = True
                i += 1
                dict_mols[i] = [line]
                path_mol = Path("toppar") / f"molecule_{i}.itp"
                fout.write(f'#include "{str(path_mol)}"\n')
                continue  # skip to next line
            if "[ defaults ]" in line:
                inside = True
                if i != 0: i += 1  # only increment if not the first molecule
                dict_mols[i] = [line]
                path_ff = Path("toppar") / f"forcefield.itp"
                logger.debug("Writing forcefield include.")
                fout.write(f'#include "{str(path_ff)}"\n')
                continue  # skip the original line
            if "[ system ]" in line or "#include" in line:
                inside = False
            if inside:
                dict_mols[i].append(line)
                continue  # skip to next line
            # write all other lines
            fout.write(line)
    # Replace the original topology with the updated file
    shutil.move(temp_file, top)
 
    for key,value in dict_mols.items():
        if "[ defaults ]" in value[0]:
            itp = path_toppar / "forcefield.itp"
        else:
            itp = path_toppar / f"molecule_{key}.itp"
        with open(itp, "w") as f:
            f.writelines(value)


    ##### CHECK OUTPUT FILES #####
    logger.info(f"The following .itp files have been created in {path_toppar} and included in {top}:")
    itps = list(path_toppar.glob(f"*.itp"))
    for itp_file in itps:
        logger.info(f"generated {itp_file}")
    
    if len(itps) == expected_itps:
        logger.info(f"All expected .itp files have been generated successfully.")
    else:
        logger.warning(f"Expected {expected_itps} .itp files, but found {len(itps)}.")
        
    # Create checkpoint file
    checkpoint = path_toppar / f"check_03_extractitp_{int(time.time())}.done"
    checkpoint.parent.mkdir(exist_ok=True, parents=True)
    with checkpoint.open('w') as f:
        f.write(".itp extraction completed successfully.\n")

if __name__ == "__main__":
    main()


    