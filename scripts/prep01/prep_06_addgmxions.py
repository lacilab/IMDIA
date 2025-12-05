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
from scripts.utils.helpers import setup_logger

def main():
    parser = argparse.ArgumentParser(description="Add GROMACS ions to the system to neutralize it and reach target concentration.")
    parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
    parser.add_argument("--gmxdir", '-g', dest='gmxdir', type=str, required=True,
                        help='Path to the GROMACS directory (containing gromacs.gro and topol.top)')
    parser.add_argument("--ions", "-i", dest='ions', type=str, default="scripts/prep01/ions.mdp", help="Path to the GROMACS parameter file to add ions.")

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    ##### Initialize logger
    logger = setup_logger("prep", args.logfile)

    # Check if checkpoint file from previous run exists
    if glob.glob(f"{args.gmxdir}/toppar/check_06_addgmxions_*.done"):
            logger.warning("Checkpoint file found.")
            # check if ions.gro exists
            if os.path.isfile(f"{args.gmxdir}/ions.gro"):
                logger.info(f"Updated topology file found at {args.gmxdir}/ions.gro.")
                logger.warning("Skipping GROMACS ions addition.")
                sys.exit(0)
            else:
                logger.warning(f"Checkpoint file found but ions.gro not found at {args.gmxdir}. The previous run might have been incomplete. Proceeding with GROMACS ions addition.")

    ##### CHECK INPUT FILES #####
    if args.ions and (not os.path.isfile(args.ions) or not args.ions.endswith(".mdp")):
        logger.error(f"Provided parameter file '{args.ions}' does not exist or is not an .mdp file.")
        sys.exit(1)

    ##### COPY INPUT FILES #####
    shutil.copy(args.ions, f"{args.gmxdir}/ions.mdp")
    # Check the ions.mdp is in the target folder
    if os.path.isfile(f"{args.gmxdir}/ions.mdp"):
        logger.info(f"Copied {args.ions} to {args.gmxdir}/ions.mdp.")
    else:
        logger.error(f"Failed to copy {args.ions} to {args.gmxdir}/ions.mdp.")
        sys.exit(1)

    ##### ADD IONS WITH GROMACS #####
    logger.debug("Adding ions to neutralize the system and reach the target concentration.")
    cmd = ["gmx", "grompp", 
           "-f", f"{args.gmxdir}/ions.mdp", 
           "-c", f"{args.gmxdir}/updated_gromacs.gro", 
           "-p", f"{args.gmxdir}/topol.top", 
           "-o", f"{args.gmxdir}/ions.tpr", 
           "-quiet", "-maxwarn", "1"]
    
    process = subprocess.run(cmd, check=True, capture_output=True, text=True)

    if process.returncode != 0:
        logger.error(f"gmx grompp error: {process.stderr}")
    else:
        logger.info(f"gmx grompp output: {process.stdout}")

    # There will be a warning because the system is not neutralized.
    cmd = ["gmx", "genion", 
           "-s", f"{args.gmxdir}/ions.tpr", 
           "-o", f"{args.gmxdir}/ions.gro", 
           "-p", f"{args.gmxdir}/topol.top", 
           "-pname", "Na+", 
           "-nname", "Cl-", 
           "-neutral", "-conc", "0.15", 
           "-quiet"]

    process = subprocess.run(
            cmd,
            input="Water\nquit\n",   # interactive stdin
            text=True,
            capture_output=True,
            check=True
        )
    
    if process.returncode != 0:
        logger.error(f"gmx genion error: {process.stderr}")
    else:
        logger.debug(f"gmx genion output: {process.stdout}")

    ##### CHECKS #####
    # Print last log line as check
    #with open(logfile, "rb") as f:
    #    lines = f.readlines()[::-1]
    #    if len(lines) == 0:
    #        logger.warning(f"Log file {logfile} is empty. There might have been an issue with the genion command.")
    #    else:
    #        for line in lines:
    #            line = line.strip()
    #            if line:
    #                logger.debug(line)
    #                break

    # check if topol.top contains more than 1 Na and Cl ions
    with open(f"{args.gmxdir}/topol.top", "r") as f:
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
    with fileinput.FileInput(f"{args.gmxdir}/topol.top", inplace=True, backup='.bak') as f:
        for line in f:
            replacements = {"Na+": "Na", "Cl-": "Cl"}
            for old, new in replacements.items():
                if old in line:
                    line = line.replace(old, new)
            print(line, end='')

    # Move gmx log file to gmxdir
    if os.path.isfile("mdout.mdp"):
        shutil.move("mdout.mdp", f"{args.gmxdir}/mdout.mdp")

    # Create checkpoint file
    checkpoint_filename = f"{args.gmxdir}/check_06_addgmxirons_{int(time.time())}.done"
    with open(checkpoint_filename, 'w') as f:
        f.write("GROMACS ions addition completed successfully.\n")

if __name__ == "__main__":
    main()