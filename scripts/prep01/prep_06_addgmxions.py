#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import fileinput
import os
import re
import subprocess
import sys
import glob
import shutil
import time
from scripts.utils.helpers import count_lines, setup_logger
from pathlib import Path

def main(logfile: Path, gmxdir: Path, ions: Path = Path("scripts/prep01/ions.mdp")):
    """
    Add GROMACS ions to the system to neutralize it and reach target concentration.

    Arguments:
        logfile : Path to the log file.
        gmxdir  : Path to the GROMACS directory (containing gromacs.gro and topol.top)
        ions    : Path to the GROMACS parameter file to add ions (default: scripts/prep01/ions.mdp)
    """

    ##### Initialize logger
    logger = setup_logger("prep", logfile)

    path_toppar = gmxdir / "toppar"

    # Check if checkpoint file from previous run exists
    checkpoint_files = list(gmxdir.glob("check_06_addgmxions_*.done"))
    if len(checkpoint_files) > 0:
        logger.warning("Checkpoint file found.")
        # check if ions.gro exists
        if os.path.isfile(gmxdir / "ions.gro"):
            logger.info(f"Found updated topology file {str(gmxdir / 'ions.gro')}.")
            logger.warning("Skipping GROMACS ions addition.")
            return
        else:
            logger.warning(f"Checkpoint file found but ions.gro not found at {str(gmxdir)}. The previous run might have been incomplete. Proceeding with GROMACS ions addition.")

    ##### CHECK INPUT FILES #####
    if ions and (not ions.is_file() or not ions.suffix == ".mdp"):
        logger.error(f"Provided parameter file '{ions}' does not exist or is not an .mdp file.")
        sys.exit(1)

    ##### COPY INPUT FILES #####
    shutil.copy(ions, gmxdir / "ions.mdp")
    # Check the ions.mdp is in the target folder
    if (gmxdir / "ions.mdp").is_file():
        logger.info(f"Copied {str(ions)} to {str(gmxdir / 'ions.mdp')}.")
    else:
        logger.error(f"Failed to copy {str(ions)} to {str(gmxdir / 'ions.mdp')}.")
        sys.exit(1)

    ##### ADD IONS WITH GROMACS #####
    logger.debug("Adding ions to neutralize the system and reach the target concentration.")
    cmd = ["gmx", "grompp", 
           "-f", f"{(gmxdir / 'ions.mdp').as_posix()}", 
           "-c", f"{(gmxdir / 'updated_gromacs.gro').as_posix()}", 
           "-p", f"{(gmxdir / 'topol.top').as_posix()}", 
           "-o", f"{(gmxdir / 'ions.tpr').as_posix()}", 
           "-po", f"{(gmxdir / 'ions_out.mdp').as_posix()}", 
           "-quiet", "-maxwarn", "1"]
    
    process = subprocess.run(cmd, capture_output=True, text=True)

    if process.returncode != 0:
        logger.error(f"gmx grompp error: {process.stderr}")
    else:
        logger.info(f"gmx grompp output: {process.stdout}")

    # # Change the name of the solvent from WAT to Water in .top so it can be selected as a group.
    replacements = {
        "WAT": "Water"
    }

    # Change ion names to gromacs standards on all .itp files in $gmxdir/toppar
    topfiles = glob.glob(f"{gmxdir}/toppar/*.itp") + glob.glob(f"{gmxdir}/topol.top")
    for file in topfiles:
        logger.debug(f"Updating solvent names in {file}...")
        num_WAT, num_water = count_lines(file, "WAT", "Water")
        logger.debug(f"Before: WAT: {num_WAT}, Water: {num_water}")
        with fileinput.FileInput(file, inplace=True) as f:
            for line in f:
                for old, new in replacements.items():
                    if old in line:
                        line = line.replace(old, new)
                print(line, end="")
        num_WAT, num_water = count_lines(file, "WAT", "Water")
        logger.debug(f"After: WAT: {num_WAT}, Water: {num_water}")
        if num_WAT > 0:
            logger.warning(f"Substitution incomplete in {file}. This could result in clashing nomenclature later.")

    # There will be a warning because the system is not neutralized.
    cmd = ["gmx", "genion", 
           "-s", f"{(gmxdir / 'ions.tpr').as_posix()}", 
           "-o", f"{(gmxdir / 'ions.gro').as_posix()}", 
           "-p", f"{(gmxdir / 'topol.top').as_posix()}", 
           "-pname", "Na+", 
           "-nname", "Cl-", 
           "-neutral", "-conc", "0.15", 
           "-quiet"]

    process = subprocess.run(
            cmd,
            input="Water\nquit\n",   # interactive stdin
            text=True,
            capture_output=True,
        )

    if process.returncode != 0:
        logger.error(f"gmx genion error: {process.stderr}")
    else:
        logger.debug(f"gmx genion output: {process.stdout}")
        
    ##### CHECKS #####
    
    # check if topol.top contains more than 1 Na and Cl ions
    with open(f"{gmxdir}/topol.top", "r") as f:
        num_Na = None
        num_Cl = None
        pattern_Na = re.compile(r"^Na\+?\b")
        pattern_Cl = re.compile(r"^Cl-?\b")
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            if pattern_Na.match(parts[0]):
                num_Na = parts[1]
            elif pattern_Cl.match(parts[0]):
                num_Cl = parts[1]

        if num_Na and int(num_Na) >= 1:
            logger.info(f"Ions added successfully. Na: {num_Na}")
        else:
            logger.warning(f"The number of Na ions added to the system is smaller or equal to 1. The system is likely not neutralized. Na: {num_Na}")

        if num_Cl and int(num_Cl) >= 1:
            logger.info(f"Ions added successfully. Cl: {num_Cl}")
        else:
            logger.warning(f"The number of Cl ions added to the system is smaller or equal to 1. The system is likely not neutralized. Cl: {num_Cl}")

    # Unify nomenclature to avoid clashes
    for file in topfiles:
        with fileinput.FileInput(file, inplace=True) as f:
            for line in f:
                replacements = {"Na+": "Na", "Cl-": "Cl", "Water": "WAT"}
                for old, new in replacements.items():
                    if old in line:
                        line = line.replace(old, new)
                print(line, end='')

    # Create checkpoint file
    checkpoint_filename = gmxdir / f"check_06_addgmxions_{int(time.time())}.done"
    with open(checkpoint_filename, 'w') as f:
        f.write("GROMACS ions addition completed successfully.\n")

if __name__ == "__main__":
    main()