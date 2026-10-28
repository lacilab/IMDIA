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
        logger.debug(f"Checkpoint file(s) found: {', '.join([str(f) for f in checkpoint_files])}")
        # Check if expected output files exist
        missing = [p for p in expected_files if not p.is_file()]
        if not missing:
            logger.warning("The expected output files are already present. Skipping Amber to GROMACS conversion.")
            return

    # If checkpoint files are not found, proceed with the conversion
    amber = pmd.load_file(str(prmtop), str(inpcrd))

    # One CA atom type per CMAP grid, so GROMACS finds each residue's grid by atom type
    _retype_cmap_atoms(amber, logger)

    # output
    amber.save(str(top), overwrite=True, format='gromacs')
    amber.save(str(gro), overwrite=True, format='gro')

    logger.debug("Checking generated files.")

    missing = [p for p in expected_files if not p.is_file()]
    for file in expected_files:
        if file not in missing:
            logger.debug(f"generated {file}")
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

def _retype_cmap_atoms(amber, logger) -> dict[str, list[str]]:
    """
    Give the central (CA) atom of each CMAP term an atom type specific to its CMAP grid.

    GROMACS looks up CMAP grids in [ cmaptypes ] by the five atom types (C N CA C N),
    while Amber force fields assign grids per residue. Each distinct pair (grid, original
    CA atom type) gets a copy of the CA atom type named XC0, XC1, ... (in order of first
    appearance, skipping names already in use), with the same mass, charge and
    Lennard-Jones parameters. The grid of each term is read from the prmtop.

    Returns {new type name: sorted residue names using it}.
    """
    if not amber.cmaps:
        logger.debug("No CMAP terms in the topology. CA atom types are unchanged.")
        return {}

    # Atom type names already present, so new names never collide with them.
    used_names = {atom.type for atom in amber.atoms}
    new_types = {}      # (grid index, original CA type name) -> new AtomType
    residues = {}       # new type name -> residue names
    grids = {}          # new type name -> grid index
    counter = 0

    for cmap in amber.cmaps:
        ca = cmap.atom3
        atom_type = ca.atom_type
        key = (amber.cmap_types.index(cmap.type), atom_type.name)

        # First term with this grid and CA type: create a copy of the CA type.
        if key not in new_types:
            while f"XC{counter}" in used_names:
                counter += 1
            name = f"XC{counter}"
            used_names.add(name)
            new_type = pmd.AtomType(name, atom_type.number, atom_type.mass,
                                    atom_type.atomic_number, name, atom_type.charge)
            # Same Lennard-Jones parameters as the original CA type.
            new_type.epsilon = atom_type.epsilon
            new_type.rmin = atom_type.rmin
            new_types[key] = new_type
            residues[name] = set()
            grids[name] = key[0]

        # Assign the CA atom to the type of its grid.
        new_type = new_types[key]
        ca.atom_type = new_type
        ca.type = new_type.name
        residues[new_type.name].add(ca.residue.name)

    # Report which residues share each new type.
    summary = {name: sorted(names) for name, names in residues.items()}
    for name, names in summary.items():
        logger.debug(f"{len(summary)} CMAP atom types created for {len(amber.cmaps)} CMAP terms.")
        logger.debug(f"CMAP atom type {name}: {', '.join(names)} (grid {grids[name]})")
    return summary


if __name__ == "__main__":
    main()