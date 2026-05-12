#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import glob
import time
import parmed as pmd
import argparse
from scripts.utils.helpers import setup_logger
from pathlib import Path

def main(logfile: Path, prmtop: Path, inpcrd: Path, top: Path, gro: Path):

    """
    Convert Amber files to GROMACS format.

    Arguments:
        logfile     : Path to the log file
        prmtop      : Path to the Prmtop file
        inpcrd      : Path to the Inpcrd file 
        top         : Path to the GROMACS topology (.top) file 
        gro         : Path to the GROMACS structure (.gro) file
    """
 
    gmxdir = top.parent

    ##### Initialize logger
    logger = setup_logger("prep", logfile)
    logger.debug(f"Starting Amber to GROMACS conversion with prmtop: {str(prmtop)}, inpcrd: {str(inpcrd)}, top: {str(top)}, gro: {str(gro)}")

    ##### Check for checkpoint files from previous runs #####
    expected_files = [ top, gro ]

    checkpoint_files = list(gmxdir.glob(f"check_02_amb2gro_*.done"))
    if len(checkpoint_files) > 0:
        logger.warning(f"Checkpoint file(s) found: {', '.join([str(f) for f in checkpoint_files])}")
        # Check if expected output files exist
        missing = [p for p in expected_files if not p.is_file()]
        if not missing:
            logger.warning("The expected output files are already present. Skipping Amber to GROMACS conversion.")
            return

    # If checkpoint files are not found, proceed with the conversion
    amber = pmd.load_file(str(prmtop), str(inpcrd))

    # Cmap correction, including PTR
    CX_dir = { 'ALA' : 'XC0', 'ARG' : 'XC1', 'ASH' : 'XC2', 'ASN' : 'XC3', 'ASP' : 'XC4', 'CYM' : 'XC5', 'CYS' : 'XC6', 
            'CYX' : 'XC5', 'GLH' : 'XC7', 'GLN' : 'XC8', 'GLU' : 'XC9', 'GLY' : 'XC10', 'HID' : 'XC5', 'HIE' : 'XC5',
            'HIP' : 'XC5', 'HYP' : 'XC5', 'ILE' : 'XC11', 'LEU' : 'XC5', 'LYN' : 'XC12', 'LYS' : 'XC12', 'MET' : 'XC5',
            'PHE' : 'XC5', 'PRO' : 'XC13', 'SER' : 'XC14', 'THR' : 'XC15', 'TRP' : 'XC5', 'TYR' : 'XC5', 'VAL' : 'XC11',
                'PTR' : 'XC5' } 

    known_types = {}
    known_types_list = []
    for cmap in amber.cmaps:  # loop through all cmap terms
        atom_type = cmap.atom3.atom_type # get atom type
        new_typename = CX_dir[cmap.atom3.residue.name] # get new CA type name (e.g. XC0)
        if new_typename in known_types_list: # if new type name in known types, directly assign atom to new type
            cmap.atom3.atom_type = known_types[cmap.type]
            cmap.atom3.type = known_types[cmap.type].name
            continue # skip

        # if new type name is not in known types, create new type
        new_type = pmd.AtomType(new_typename, atom_type.number, atom_type.mass, atom_type.atomic_number, new_typename, atom_type.charge)
        new_type.epsilon = atom_type.epsilon # copy over epsilon and rmin
        new_type.rmin = atom_type.rmin
        known_types[cmap.type] = new_type
        known_types_list.append(new_typename)

        # assgin atom to new type
        cmap.atom3.atom_type = known_types[cmap.type]
        cmap.atom3.type = known_types[cmap.type].name


    # output 
    amber.save(str(top), overwrite=True, format='gromacs')
    amber.save(str(gro), overwrite=True, format='gro')

    logger.debug("Checking generated files.")

    missing = [p for p in expected_files if not p.is_file()]
    for file in expected_files:
        if file not in missing:
            logger.info(f"generated {file}")
        else:
            logger.error(f"{file} not generated")

    if not missing:
        logger.info("Amber to GROMACS conversion completed successfully.")
        checkpoint = gmxdir / f"check_02_amb2gro_{int(time.time())}.done"
        checkpoint.parent.mkdir(exist_ok=True, parents=True)
        with checkpoint.open('w') as f:
            f.write("Amber to GROMACS conversion completed successfully.\n")
        return
    else:
        logger.error("Amber to GROMACS conversion failed. There are missing files.")
        sys.exit(1)

if __name__ == "__main__":
    main()