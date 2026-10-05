#!/usr/bin/env python
# -*- coding: utf-8 -*-
from datetime import datetime
import sys
import os
import re
import subprocess
import shutil
from pathlib import Path
from scripts.utils.helpers import setup_logger
from scripts.prep01 import prep_01_amber as amber
from scripts.prep01 import prep_02_amb2gro as amb2gro
from scripts.prep01 import prep_03_extractitp as extractitp
from scripts.prep01 import prep_04_removeamberions as removeamberions
from scripts.prep01 import prep_05_restraindnahbonds as restraindnahbonds
from scripts.prep01 import prep_06_addgmxions as addgmxions
from scripts.prep01.prep_00_protonation import normalize_config

def main(logfile: Path, datadir: Path, filename: str,
    dna : bool = False, phos : bool = False,
    range_a : str = None, range_b : str = None,
    ffs: str = None, protonation: dict = None):
    """
    Prepare a PDB structure for MD simulations.

    Steps:
        1.Amber preparation: setup forcefields, add hydrogens.
        2.Conversion from Amber to GROMACS.
        3.Extract .itp files from .top topology file.
        4.Remove ions added with Amber.
        5.Add ions with GROMACS to neutralize the system and reach target concentration.

    Arguments:
        logfile  : Path to the log file.
        datadir  : Directory containing the input PDB file.
        filename : Name of the PDB file (with or without .pdb)
        dna      : True if the system contains DNA. Adds DNA.bsc1 force field.
        phos     : True if the system contains phosphorylated residues.
        range_a  : Range of residues for DNA A (0-indexed) (default: [0, -1], full length)
        range_b  : Range of residues for DNA B (0-indexed) (default: [0, -1], full length)
        ffs      : List of force fields to use in LEaP (default: protein.ff19SB, water.opc)
        protonation : Protonation configuration (default: None)
    """

    ##### INITIALIZE LOG
    logger = setup_logger("prep", logfile)

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

    protonation = normalize_config(protonation)

    ##### CHECK KEY SOFTWARE DEPENDENCIES ######
    # Check GROMACS is ready and print version
    if shutil.which("gmx") is None:
        logger.error("GROMACS (gmx) not found. Please install GROMACS and ensure it is in your PATH.", file=sys.stderr)
        sys.exit(1)
    subprocess.run(["gmx", "--version"])

    # Check ambpdb is ready and print version
    if shutil.which("ambpdb") is None:
        logger.error("ERROR: ambpdb not found. Please install ambpdb and ensure it is in your PATH.", file=sys.stderr)
        sys.exit(1)
    subprocess.run(["ambpdb", "--version"])

    # Check ParmEd is ready and print version
    if shutil.which("parmed") is None:
        logger.error("ERROR: ParmEd (parmed) not found. Please install ParmEd and ensure it is in your PATH.", file=sys.stderr)
        sys.exit(1)
    subprocess.run(["parmed", "--version"])

    logger.info(f"Preparing {datadir}/{file}.pdb")

    # Check if datadir/amber already exists.
    path_amber = datadir / "amber"
    if path_amber.is_dir():
        logger.warning(f"{str(path_amber)} already exists.")
    else:
        logger.debug(f"{str(datadir)}/amber does not exist. Creating new directory.")
        path_amber.mkdir(exist_ok=True)

    # Remove lines starting with CONECT, which are added by MOE and cause issues in AmberTools 
    path_noconect = Path(f"{datadir}/amber/{file}_noconect.pdb")
    if path_noconect.is_file():
        logger.warning(f"{str(path_noconect)} already exists. Skipping CONECT removal.")
    else:
        connect_before = 0
        with open(path_pdb) as f:
            connect_before = sum(1 for line in f if line.startswith("CONECT"))
        if connect_before == 0:
            logger.info(f"No CONECT lines found in {file}.pdb. No need to remove them.")
            shutil.copy(path_pdb, path_noconect)
            logger.debug(f"Copied {file}.pdb to amber/{file}_noconect.pdb without changes.")
        else:
            logger.debug(f"Removing CONECT lines from {file}.pdb. These lines can cause issues in AmberTools.")
            with open(path_pdb) as fin, open(path_noconect, "w") as fout:
                fout.writelines(line for line in fin if not line.startswith("CONECT"))
            remaining = 0
            with open(path_noconect) as f:
                remaining = sum(1 for line in f if line.startswith("CONECT"))
            logger.debug(f"CONECT lines before: {connect_before}")
            logger.debug(f"CONECT lines after:  {remaining}")
            if remaining > 0:
                logger.warning(f"{remaining} CONECT lines remain in {file}_noconect.pdb")

    if path_amber.is_dir() and path_noconect.is_file():
        logger.info(f"{str(path_amber)} directory is ready. ")
        logger.info(f"{file}_noconect.pdb is in {str(path_amber)}/")
    else:
        logger.error(f"{str(path_amber)} directory is missing or incomplete.")
        sys.exit(1)

    ##### PREPARATION ######

    ### 1. AmberTools preparation

    amber.main(logfile=logfile, pdbpath=path_noconect, dna=dna, phos=phos, ffs=ffs, protonation=protonation)

    logger.debug(f"The current directory is {os.getcwd()}")

    # Check that datadir/gromacs_common does not exist. If it exists, show a message and backup the old one (change the name)
    gmxdir = datadir / "gromacs_common"
    if gmxdir.is_dir():
        logger.warning(f"{str(gmxdir)} already exists.")
    else:
        gmxdir.mkdir(exist_ok=True)
        if gmxdir.is_dir():
            logger.info(f"GROMACS directory created at {str(gmxdir)}.")
        else:
            logger.error(f"Failed to create GROMACS directory at {str(gmxdir)}.")
            sys.exit(1)

    ### 2. Conversion from Amber to GROMACS with parmed

    path_prmtop = datadir / "amber" / f"{file}_noconect_protonated_solvated.prmtop"
    path_inpcrd = datadir / "amber" / f"{file}_noconect_protonated_solvated.inpcrd"
    path_top = gmxdir / "topol.top"
    path_gro = gmxdir / "gromacs.gro"
    path_toppar = gmxdir / "toppar"

    amb2gro.main(logfile=logfile, prmtop=path_prmtop,
        inpcrd=path_inpcrd,
        top=path_top,
        gro=path_gro)

    ### 3. Extract .itps from .top file
    extractitp.main(logfile=logfile, top=path_top, dna=dna)

    ### 4. Remove the 2 ions added with AmberTools from .gro files
    removeamberions.main(logfile=logfile, gmxdir=gmxdir)

    ### 5. Add ions with GROMACS (neutralize and reach target concentration)
    addgmxions.main(logfile=logfile, gmxdir=gmxdir)

    ### 6. Add DNA H-bonds restraints if DNA is present
    if dna:
        path_top = gmxdir / "topol.top"
        restraindnahbonds.main(logfile=logfile, gmxdir=gmxdir, topol=path_top, range_a=range_a, range_b=range_b)
    else:
        logger.debug("The complex does not include DNA. Skipping DNA H-bond restraints. If this is incorrect, rerun with dna: True in params.yml.")

    logger.info("Preparation completed successfully.")


if __name__ == "__main__":
    main()