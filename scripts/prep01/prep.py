#!/usr/bin/env python
# -*- coding: utf-8 -*-
from datetime import datetime
import sys
import os
import glob
import re
import subprocess
import shutil
import argparse
from scripts.utils.helpers import setup_logger

parser = argparse.ArgumentParser(description="Prepare a PDB structure for MD simulations.")
parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
parser.add_argument('--datadir', '-g', dest='datadir', default='.', help='Directory that contains the input PDB file.')
parser.add_argument('--filename', '-f', dest='filename', default=None, help='Name of the PDB file (with or without .pdb).')
parser.add_argument('--dna', dest='with_dna', action='store_true', help='Use if the system contains DNA. Adds DNA.bsc1 force field.')
parser.add_argument('--phos', dest='with_phosphate', action='store_true', help='Use if the system contains phosphorylated residues.')
parser.add_argument('--range_a', dest='range_a', type=str, default=None, help='Range of residues to consider for DNA strand A (e.g., 1,20). Default is full length.')
parser.add_argument('--range_b', dest='range_b', type=str, default=None, help='Range of residues to consider for DNA strand B (e.g., 1,20). Default is full length.')
parser.add_argument('--ffs', dest='ffs', type=str, default=[], help='Space-separated list of force fields to use in LEaP. Default: protein.ff19SB water.opc.')

# print help if required arguments are not provided
if len(sys.argv) < 3:
    parser.print_help()
    sys.exit(0)

args = parser.parse_args()

##### INITIALIZE LOG
logger = setup_logger("prep", args.logfile)

###### CHECK INPUTS ######
# Normalize filename: accept with or without .pdb, store basename without extension
file = os.path.basename(args.filename)
file, ext = os.path.splitext(file)

# check datadir exists
if not os.path.isdir(args.datadir):
    logger.error(f"datadir '{args.datadir}' does not exist or is not a directory.")
    sys.exit(1)
# remove trailing / from datadir
datadir = args.datadir.rstrip("/")

# check input PDB file exists
if not os.path.isfile(f"{datadir}/{file}.pdb"):
    logger.error(f"Input PDB file '{datadir}/{file}.pdb' not found.")
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
if os.path.isdir(f"{datadir}/amber"):
    logger.warning(f"{datadir}/amber already exists.")
else:
    print(f"{datadir}/amber does not exist. Creating new directory.")
    os.makedirs(f"{datadir}/amber")


# Remove lines starting with CONECT, which are added by MOE and cause issues in AmberTools 
if os.path.isfile(f"{datadir}/amber/{file}_noconect.pdb"):
    logger.warning(f"{datadir}/amber/{file}_noconect.pdb already exists. Skipping CONECT removal.")
else:
    connect_before = 0
    with open(f"{datadir}/{file}.pdb") as f:
        connect_before = sum(1 for line in f if line.startswith("CONECT"))
    if connect_before == 0:
        logger.info(f"No CONECT lines found in {file}.pdb. No need to remove them.")
        shutil.copy(f"{datadir}/{file}.pdb", f"{datadir}/amber/{file}_noconect.pdb")
        print(f"Copied {file}.pdb to amber/{file}_noconect.pdb without changes.")
    else:
        print(f"Removing CONECT lines from {file}.pdb. These lines can cause issues in AmberTools.")
        with open(f"{datadir}/{file}.pdb") as fin, open(f"{datadir}/amber/{file}_noconect.pdb", "w") as fout:
            fout.writelines(line for line in fin if not line.startswith("CONECT"))
        remaining = 0
        with open(f"{datadir}/amber/{file}_noconect.pdb") as f:
            remaining = sum(1 for line in f if line.startswith("CONECT"))
        print(f"CONECT lines before: {connect_before}")
        print(f"CONECT lines after:  {remaining}")
        if remaining > 0:
            logger.warning(f"{remaining} CONECT lines remain in {file}_noconect.pdb")

if os.path.isdir(f"{datadir}/amber") and os.path.isfile(f"{datadir}/amber/{file}_noconect.pdb"):
    logger.info(f"{datadir}/amber directory is ready. ")
    logger.info(f"{file}_noconect.pdb is in {datadir}/amber/")
else:
    logger.error(f"{datadir}/amber directory is missing or incomplete.")
    sys.exit(1)


##### PREPARATION ######

### 1. AmberTools preparation

dna_flag = "--dna" if args.with_dna else ""
phos_flag = "--phos" if args.with_phosphate else ""
ffs_flag = "--ffs " if args.ffs else ""
ffs_val = ','.join([str(x) for x in args.ffs]) if args.ffs else ""

cmd = [sys.executable, "-m", "scripts.prep01.prep_01_amber", 
        "--logfile", f"{args.logfile}",
       "--path-pdb", f"{datadir}/amber/{file}_noconect.pdb",
       dna_flag,
       phos_flag,
       ffs_flag, ffs_val,
    ]
cmd = [x for x in cmd if x.strip()]
subprocess.run(cmd, check=True)

logger.debug(f"The current directory is {os.getcwd()}")

# Check that datadir/gromacs_common does not exist. If it exists, show a message and backup the old one (change the name)
gmxdir = f"{datadir}/gromacs_common"
if os.path.isdir(gmxdir):
    logger.warning(f"{gmxdir} already exists.")
else:
    os.makedirs(gmxdir)
    if os.path.isdir(f"{datadir}/gromacs_common"):
        logger.info(f"GROMACS directory created at {gmxdir}.")
    else:
        logger.error(f"Failed to create GROMACS directory at {gmxdir}.")
        sys.exit(1)

### 2. Conversion from Amber to GROMACS with parmed

cmd = [sys.executable, "-m", "scripts.prep01.prep_02_amb2gro", 
        "--logfile", f"{args.logfile}",
       "-p", f"{datadir}/amber/{file}_noconect_protonated_solvated.prmtop",
        "-c", f"{datadir}/amber/{file}_noconect_protonated_solvated.inpcrd",
        "-t", f"{gmxdir}/topol.top",
        "-g", f"{gmxdir}/gromacs.gro"]

subprocess.run(cmd, check=True)

### 3. Extract .itps from .top file
cmd = [sys.executable, "-m", "scripts.prep01.prep_03_extractitp",  "--logfile", f"{args.logfile}", "--top", f"{gmxdir}/topol.top", dna_flag]
subprocess.run(cmd, check=True)

### 4. Remove the 2 ions added with AmberTools from .gro files
cmd = [sys.executable, "-m", "scripts.prep01.prep_04_removeamberions", "--logfile", f"{args.logfile}","--gmxdir", f"{gmxdir}"]
subprocess.run(cmd, check=True)

### 5. Add ions with GROMACS (neutralize and reach target concentration)
cmd = [sys.executable, "-m", "scripts.prep01.prep_06_addgmxions", "--logfile", f"{args.logfile}", "--gmxdir", f"{gmxdir}"]
subprocess.run(cmd, check=True)

### 6. Add DNA H-bonds restraints if DNA is present
if args.with_dna:
    dna_files = []
    for file in glob.glob(f"{gmxdir}/toppar/"+"*.itp"):
        if re.search(r'DC|DA|DG|DT', open(file).read()):
            print(f"Found DNA base in {file}")
            dna_files.append(file)
    if len(dna_files) != 2:
        logger.error(f"--dna flag provided but found {len(dna_files)} files containing DNA residues. Expected 2 files.")
        sys.exit(1)
    print("The complex includes DNA. Adding DNA H-bond restraints.")

    flag_range_a = ["--range_a", args.range_a] if args.range_a else []
    flag_range_b = ["--range_b", args.range_b] if args.range_b else []

    cmd = [sys.executable, "-m", "scripts.prep01.prep_05_restraindnahbonds", "--logfile", f"{args.logfile}", "--gmxdir", gmxdir,
        "--topol", f"{gmxdir}/topol.top",
        "--dna_a", dna_files[0],
        "--dna_b", dna_files[1]] + flag_range_a + flag_range_b
    subprocess.run(cmd, check=True, text = True)
else:
    print("The complex does not include DNA. Skipping DNA H-bond restraints. If this is incorrect, rerun with the --dna flag.")

logger.info("Preparation completed successfully.")