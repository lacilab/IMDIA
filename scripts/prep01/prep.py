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

# parser = argparse.ArgumentParser(description="Prepare a PDB structure for MD simulations.")
# parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
# parser.add_argument('--datadir', '-g', dest='datadir', default='.', help='Directory that contains the input PDB file.')
# parser.add_argument('--filename', '-f', dest='filename', default=None, help='Name of the PDB file (with or without .pdb).')
# parser.add_argument('--dna', dest='with_dna', action='store_true', help='Use if the system contains DNA. Adds DNA.bsc1 force field.')
# parser.add_argument('--phos', dest='with_phosphate', action='store_true', help='Use if the system contains phosphorylated residues.')
# parser.add_argument('--range_a', dest='range_a', type=str, default=None, help='Range of residues to consider for DNA strand A (e.g., 1,20). Default is full length.')
# parser.add_argument('--range_b', dest='range_b', type=str, default=None, help='Range of residues to consider for DNA strand B (e.g., 1,20). Default is full length.')
# parser.add_argument('--ffs', dest='ffs', type=str, default=[], help='Space-separated list of force fields to use in LEaP. Default: protein.ff19SB water.opc.')

def main(logfile: Path, datadir: Path, filename: str,
    dna : bool = False, phos : bool = False,
    range_a : str = None, range_b : str = None,
    ffs : str = None):

    ##### INITIALIZE LOG
    logger = setup_logger("prep", logfile)

    ###### CHECK INPUTS ######
    # Normalize filename: accept with or without .pdb, store basename without extension
    file = Path(filename).stem

    # check datadir exists
    if not datadir.is_dir():
        logger.error(f"datadir '{datadir}' does not exist or is not a directory.")
        sys.exit(1)
    # remove trailing / from datadir
    # datadir = datadir.rstrip("/")
    
    # check input PDB file exists
    path_pdb = datadir / f"{file}.pdb"
    if not path_pdb.is_file():
        logger.error(f"Input PDB file '{str(path_pdb)}' not found.")
        sys.exit(1)

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

    # Check for checkpoint files
    #if glob.glob(f"{datadir}/check_01_prep_*.done"):
    #    logger.warning("Checkpoint file found. Skipping preparation.")
    #    print(f"To rerun these steps, please remove the checkpoint file(s): {datadir}/check_01_prep_*.done")
    #    sys.exit(0)

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

    amber.main(logfile=logfile, pdbpath=path_noconect, dna=dna, phos=phos, ffs=ffs)

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
    # cmd = [sys.executable, "-m", "scripts.prep01.prep_04_removeamberions", "--logfile", f"{logfile}","--gmxdir", f"{gmxdir}"]
    # subprocess.run(cmd, check=True)

    removeamberions.main(logfile=logfile, gmxdir=gmxdir)

    ### 5. Add ions with GROMACS (neutralize and reach target concentration)
    # cmd = [sys.executable, "-m", "scripts.prep01.prep_06_addgmxions", "--logfile", f"{logfile}", "--gmxdir", f"{gmxdir}"]
    # subprocess.run(cmd, check=True)

    addgmxions.main(logfile=logfile, gmxdir=gmxdir)

    ### 6. Add DNA H-bonds restraints if DNA is present
    if dna:
        dna_files = []
        for file in path_toppar.glob("*.itp"):
            if re.search(r'DC|DA|DG|DT', open(file).read()):
                logger.debug(f"Found DNA base in {file}")
                dna_files.append(Path(file))
        if len(dna_files) != 2:
            logger.error(f"--dna flag provided but found {len(dna_files)} files containing DNA residues. Expected 2 files.")
            sys.exit(1)
        logger.info("The complex includes DNA. Adding DNA H-bond restraints.")

        path_top = gmxdir / "topol.top"
        restraindnahbonds.main(logfile=logfile, gmxdir=gmxdir, topol=path_top,
            dna_a=dna_files[0], dna_b=dna_files[1], range_a=range_a, range_b=range_b)
    else:
        print("The complex does not include DNA. Skipping DNA H-bond restraints. If this is incorrect, rerun with the --dna flag.")

    logger.info("Preparation completed successfully.")


if __name__ == "__main__":
    main()