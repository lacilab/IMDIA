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


def main():
    parser = argparse.ArgumentParser(description="Generate .itp files from topology file.")
    parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
    parser.add_argument("--top", '-p', dest='top', type=str, default='topol.top',
                        help='Path to the GROMACS topology file (e.g. topol.top)')
    parser.add_argument('--dna', dest='with_dna', action='store_true', help='Use if the system contains DNA.')

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    gmxdir = os.path.dirname(args.top)
    expected_itps = 8 if args.with_dna else 6

    ##### Initialize logger
    logger = setup_logger("prep", args.logfile)

    # Check for checkpoint files
    if glob.glob(f"{gmxdir}/toppar/check_03_extractitp_*.done"):
        logger.warning("Checkpoint file found. Checking .itp files.")
        itps = glob.glob(f"{gmxdir}/toppar/*.itp")
        if len(itps) != expected_itps:
            logger.warning(f"Expected {expected_itps} .itp files, but found {len(itps)}. The previous run might have been incomplete. Proceeding with .itp extraction.")
        else:
            logger.info(f"The number of .itp files found matches the expected count of ({expected_itps}).")
            logger.warning("Skipping .itp extraction.")
            sys.exit(0)

    # Check that gmxdir/toppar does not exist
    if not os.path.isdir(f"{gmxdir}/toppar"):
        os.makedirs(f"{gmxdir}/toppar")
        if os.path.isdir(f"{gmxdir}/toppar"):
            logger.info(f"{gmxdir}/toppar directory created at {gmxdir}.")
        else:
            logger.error(f"{gmxdir}/toppar directory could not be created at {gmxdir}.")
            sys.exit(1)
    else:
        logger.warning(f"{gmxdir}/toppar already exists")

    # Extract molecule itps
    pattern = re.compile(r"\[\s*moleculetype\s*\]")
    with open(args.top, "r") as f:
        n = sum(1 for line in f if pattern.search(line))
    logger.debug(f"Found {n} molecule types in {args.top}.")
    i = 0
    temp_file = args.top + ".tmp"
    inside = False
    dict_mols = {}
    with open(args.top, 'r') as fin, open(temp_file, "w") as fout:
        for line in fin:
            if "[ moleculetype ]" in line:
                inside = True
                i += 1
                dict_mols[i] = [line]
                fout.write(f'#include "toppar/molecule_{i}.itp"\n')
                continue  # skip to next line
            if "[ defaults ]" in line:
                inside = True
                if i != 0: i += 1  # only increment if not the first molecule
                dict_mols[i] = [line]
                logger.debug("Writing forcefield include.")
                fout.write(f'#include "toppar/forcefield.itp"\n')
                continue  # skip the original line
            if "[ system ]" in line or "#include" in line:
                inside = False
            if inside:
                dict_mols[i].append(line)
                continue  # skip to next line
            # write all other lines
            fout.write(line)
    # Replace the original topology with the updated file
    shutil.move(temp_file, args.top)
 
    for key,value in dict_mols.items():
        if "[ defaults ]" in value[0]:
            itp = f"{gmxdir}/toppar/forcefield.itp"
        else:
            itp = f"{gmxdir}/toppar/molecule_{key}.itp"
        with open(itp, "w") as f:
            f.writelines(value)


    ##### CHECK OUTPUT FILES #####
    logger.info(f"The following .itp files have been created in {gmxdir}/toppar/ and included in {args.top}:")
    for itp_file in glob.glob(f"{gmxdir}/toppar/*.itp"):
        logger.info(f"generated {itp_file}")
    
    # Check that the expected number of itp files have been created
    itps = glob.glob(f"{gmxdir}/toppar/*.itp")
    if len(itps) == expected_itps:
        logger.info(f"All expected .itp files have been generated successfully.")
        # Create checkpoint file
        checkpoint_filename = f"{gmxdir}/toppar/check_03_extractitp_{int(time.time())}.done"
        with open(checkpoint_filename, 'w') as f:
            f.write(".itp extraction completed successfully.\n")
    else:
        logger.warning(f"Expected {expected_itps} .itp files, but found {len(itps)}.")
        

if __name__ == "__main__":
    main()


    