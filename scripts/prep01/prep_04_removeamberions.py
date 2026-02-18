#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import fileinput
import os
import re
import sys
import glob
import time
from pathlib import Path
from scripts.utils.helpers import count_lines, setup_logger


def main(logfile: Path, gmxdir: Path):
    # parser = argparse.ArgumentParser(description="Remove Amber ions from the GROMACS topology files.")
    # parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
    # parser.add_argument("--gmxdir", '-g', dest='gmxdir', type=str, required=True,
    #                     help='Path to the GROMACS directory (containing gromacs.gro and topol.top)')
    # )

    ##### Initialize logger
    logger = setup_logger("prep", logfile)

    path_toppar = gmxdir / "toppar"
    gro = gmxdir / "gromacs.gro"
    updated_gro = gmxdir / "updated_gromacs.gro"
    top = gmxdir / "topol.top"

    ##### CHECK INPUT FILES #####
    if not gmxdir.is_dir():
        logger.error(f"GROMACS directory '{str(gmxdir)}' does not exist or is not a directory.")
        sys.exit(1)
    if not gro.is_file():
        logger.error( f"gromacs.gro file not found in '{str(gmxdir)}'.")
        sys.exit(1)
    if not top.is_file():
        logger.error( f"topol.top file not found in '{str(gmxdir)}'.")
        sys.exit(1)
    if not path_toppar.is_dir():
        logger.error( f"toppar directory not found in '{str(gmxdir)}'.")
        sys.exit(1)
    # check number of .itp files in toppar is greater or equal to 2
    itp_files = list(path_toppar.glob(f"*.itp"))
    if len(itp_files) < 2:
        logger.error( f"No .itp files found in '{str(path_toppar)}'.")
        sys.exit(1)

    # Check if checkpoint file from previous run exists
    if len(list(path_toppar.glob(f"check_04_removeamberions_*.done"))) > 0:
        logger.warning("Checkpoint file found.")
         # check if updated_gromacs.gro exists
        path_updated_gro = gmxdir / "updated_gromacs.gro"
        if os.path.isfile(path_updated_gro):
            logger.info(f"Updated topology file found at {str(path_updated_gro)}.")
            logger.warning("Skipping Amber ions removal.")
            return
        else:
            logger.warning(f"Checkpoint file found but updated_gromacs.gro not found at {str(gmxdir)}. The previous run might have been incomplete. Proceeding with Amber ions removal.")



    ##### HELPERS ######
    # Helper function to wrap molecule and atom numbers
    def wrap_minus_2(n):
        return n - 2 if n >= 2 else 99998 + n
    
    ##### REMOVE AMBER IONS FROM TOPOLOGY AND GRO FILES #####
    # count lines with Na and Cl in gromacs.gro (before removal)
    line_Na = None
    line_Cl = None
    with open(gro) as f:
        num_Na, num_Cl = count_lines(gro, "Na", "Cl")
        # Find first lines containing "Na" and "Cl"
        for line_no, line in enumerate(f, start=1):
            # total line count
            last_line = line_no
            # first Na
            if line_Na is None and "Na" in line:
                line_Na = line_no
            # first Cl
            if line_Cl is None and "Cl" in line:
                line_Cl = line_no
    logger.debug(f"Number of Na ions in gromacs.gro: {num_Na}")
    logger.debug(f"Number of Cl ions in gromacs.gro: {num_Cl}")

    # Warnings if num_Na or num_Cl is more than 1
    if num_Na > 1:
        logger.warning(f"Number of Na ions is expected to be one, but was {num_Na}. Later addition of ions with GROMACS may not work as expected.")
    if num_Cl > 1:
        logger.warning(f"Number of Cl ions is expected to be one, but was {num_Cl}. Later addition of ions with GROMACS may not work as expected.")

    # Check if both Na and Cl lines were found
    if line_Na is None or line_Cl is None:
        logger.error("Could not find both Na and Cl lines in gromacs.gro.")
        sys.exit(1)
    else:
        logger.info(f"Found Na at line {line_Na} and Cl at line {line_Cl} in gromacs.gro.")

    
    line1, line2 = sorted([line_Na, line_Cl])

    # Remove Na and Cl lines and save as updated_gromacs.gro
    with open(gro, "r") as fin, open(updated_gro, "w") as fout:
        for NR, line in enumerate(fin, start=1):
            # NR == 2 : update atom count
            if NR == 2:
                atom_count = int(line[:6]) - 2
                fout.write(f"{atom_count}\n")
                continue
            # Copy lines before line1
            if NR < line1:
                fout.write(line)
                continue

            # Skip removed lines
            if NR == line1 or NR == line2:
                continue

            # Last line → output untouched
            if NR == last_line:
                fout.write(line)
                continue

            # Process atom lines (NR > line2 and NR < last_line)

            # Split columns by one or more spaces
            cols = re.split(r'\s+', line.strip())

            # extraction of mol_number + mol_name
            # mol looks like "   12Na" → number+text glued
            # Find first non-digit (fast loop, avoids regex)
            s = cols[0].lstrip()
            j = 0
            while j < len(s) and s[j].isdigit():
                j += 1
            mol_number = int(s[:j])
            mol_name   = s[j:]

            # decrement
            mol_number = wrap_minus_2(mol_number)

            # Detect integer vs float atom number column
            if "." not in cols[2]:
                # integer case
                anum = wrap_minus_2(int(cols[2]))
                fout.write(
                    f"{mol_number:5d}{mol_name:<3}{cols[1]:>7}{anum:5d}"
                    f"{float(cols[3]):8.3f}{float(cols[4]):8.3f}{float(cols[5]):8.3f}\n"
                )
            else:
                # float case → atom number is embedded in atom_name
                # atom_name e.g. " H   12"
                # last 5 chars hold the number
                anum = wrap_minus_2(int(cols[1][-5:]))
                atom_str = cols[1][:-5]

                fout.write(
                    f"{mol_number:5d}{mol_name:<3}{atom_str:>7}{anum:5d}"
                    f"{float(cols[2]):8.3f}{float(cols[3]):8.3f}{float(cols[4]):8.3f}\n"
                )

    # check updated_gromacs.gro was created
    if updated_gro.is_file():
        logger.info(f"Updated gromacs.gro file created at {str(updated_gro)}.")
    else:
        logger.error(f"Failed to create updated gromacs.gro file at {str(updated_gro)}.")
        sys.exit(1)

    # check there are no Na or Cl updated_gromacs.gro
    with open(updated_gro) as f:
        num_Na_updated, num_Cl_updated = count_lines(updated_gro, "Na", "Cl")
    if num_Na_updated == 0 and num_Cl_updated == 0:
        logger.info("Successfully removed Amber ions from gromacs.gro. No Na or Cl lines found in updated_gromacs.gro.")
    else:
        logger.warning(f"Found {num_Na_updated} Na and {num_Cl_updated} Cl lines in updated_gromacs.gro, but expected 0.")

    # remove Na and Cl from topol.top
    with fileinput.FileInput(top, inplace=True) as f:
        for line in f:
            if "Na" in line or "Cl" in line:
                continue
            print(line, end="")

    replacements = {
        "Na+": "Na",
        "Cl-": "Cl",
    }

    # Change ion names to gromacs standards on all .itp files in $gmxdir/toppar
    files = list(path_toppar.glob(f"*.itp")) + [top]
    for file in files:
        logger.debug(f"Updating ion names in {file}...")
        num_Na, num_Cl, num_Nasig, num_Clsig, num_WAT, num_water = count_lines(file, "Na", "Cl", "Na+", "Cl-", "WAT", "Water")
        logger.debug(f"Before: Na+: {num_Nasig}, Cl-: {num_Clsig}, Na: {num_Na}, Cl: {num_Cl}, WAT: {num_WAT}, Water: {num_water}")
        with fileinput.FileInput(file, inplace=True) as f:
            for line in f:
                for old, new in replacements.items():
                    if old in line:
                        line = line.replace(old, new)
                print(line, end="")
        num_Na, num_Cl, num_Nasig, num_Clsig, num_WAT, num_water = count_lines(file, "Na", "Cl", "Na+", "Cl-", "WAT", "Water")
        logger.debug(f"After: Na+: {num_Nasig}, Cl-: {num_Clsig}, Na: {num_Na}, Cl: {num_Cl}, WAT: {num_WAT}, Water: {num_water}")
        if num_Nasig > 0 or num_Clsig > 0:
            logger.warning(f"Substitution incomplete in {file}. This could result in clashing nomenclature later.")

    # Create checkpoint file
    checkpoint_filename = f"{gmxdir}/check_04_removeamberions_{int(time.time())}.done"
    with open(checkpoint_filename, 'w') as f:
        f.write("Amber ions removal completed successfully.\n")


if __name__ == "__main__":
    main()
