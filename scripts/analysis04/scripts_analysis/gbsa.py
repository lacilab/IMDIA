#!/usr/bin/env python3

import sys
import os
import re
import shutil
import multiprocessing
import subprocess
from pathlib import Path
from venv import logger
from scripts.utils.helpers import setup_logger
from scripts.utils.group_index import group_index

def prep_mmpbsain(logger, infile, start, stop, dna, phos, ffs):
    """ Prepare mmpbsa.in parameter file for GBSA. """
    # get GROMACS path with which gmx. Take part before "bin"
    gmx_path = shutil.which("gmx").split("/gmx")[0]
    logger.debug(f"Detecting GROMACS path using 'which gmx' command: {gmx_path}")
    if gmx_path is None:
        logger.error("GROMACS (gmx) not found. Please install GROMACS and ensure it is in your PATH.")
        sys.exit(1)
    else:
        logger.debug(f"GROMACS found at {gmx_path}. This path will be used for preparing mmpbsa.in file.")

    if ffs is None or len(ffs) == 0:
        logger.warning("No force fields specified for GBSA calculation. Default force fields will be used.")
        ffs = ['water.opc', 'protein.ff19SB']
        if dna:
            ffs.append('dna.bsc1')
        if phos:
            ffs.append('phosaa19SB')

    # read the input file
    with open(infile, 'r+') as f:
        lines = f.readlines()
        # insert GROMACS path
        for i, line in enumerate(lines):
            if re.match(r"\s*gmx_path\s*=", line):
                lines[i] = f"  gmx_path = \"{gmx_path}\"\n"
                logger.debug(f"Inserted GROMACS path into mmpbsa.in: {gmx_path}")
            elif re.match(r"\s*forcefields\s*=", line):
                ffs_import = ",".join(f"leaprc.{ff}" for ff in ffs)
                lines[i] = f"  forcefields = \"{ffs_import}\"\n"
                logger.debug(f"Inserted force fields into mmpbsa.in: {ffs_import}")
            elif re.match(r"\s*startframe\s*=", line):
                if start is not None:
                    lines[i] = f"  startframe = {start}\n"  
                    logger.debug(f"Set startframe to {start} in mmpbsa.in")
                else:
                    lines[i] = f"  startframe = 1\n"  
                    logger.debug("No start frame provided; keeping default 1 in mmpbsa.in")
            elif re.match(r"\s*endframe\s*=", line):
                if stop is not None:
                    lines[i] = f"  endframe = {stop}\n"  
                    logger.debug(f"Set endframe to {stop} in mmpbsa.in")
                else:
                    lines[i] = f"  endframe = 9999999\n"  
                    logger.debug("No stop frame provided; keeping default 9999999 in mmpbsa.in")
        # write
        f.seek(0)
        f.writelines(lines)
        f.truncate()

def main(logfile: str = None, scriptdir: str = None, gmxdir: str = None, outdir: str = None, 
         index:str = None, tprfile: str = None,
         xtcfile: str = None, topfile: str = None, cores: int = 1,
         start: int = None, stop: int = None, 
         dna: bool = False, phos: bool = False,
         index_group1: str = None, index_group2: str = None, ffs: list = None):

    """
    Run MM/GBSA free binding energy calculations on MD trajectories.

    Arguments:
        logfile        : Path to the log file.
        scriptdir      : Directory containing mmpbsa.in.
        gmxdir         : Path to the directory with the MD results (e.g. rep1/gromacs).
        outdir         : Output directory for MM/GBSA results (e.g. rep1/analysis/gbsa).

        index          : GROMACS index file (.ndx) defining atom groups.
        tprfile        : GROMACS run input file (.tpr) containing topology and run parameters.
        xtcfile        : GROMACS trajectory file (.xtc) containing trajectory data.
        topfile        : GROMACS topology file (.top).

        cores          : Number of CPU cores used for parallel computation.

        start          : Starting frame index for trajectory analysis.
        stop           : Ending frame index for trajectory analysis.

        dna            : True if DNA is present in the system.
        phos           : True if phosphorylated residues are present.

        index_group1   : First index group for binding energy calculation.
        index_group2   : Second index group for binding energy calculation.

        ffs            : List of force fields used in the MM/GBSA setup.
    """
    # Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug("Starting GBSA analysis")

    # Create output directory if it doesn't exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        logger.info(f"Created output directory: {outdir}")
    else:
        logger.warning(f"Directory already exists: {outdir}")

    # Copy input file
    input_file_src = os.path.join(scriptdir, "mmpbsa.in")
    input_file_dst = os.path.join(outdir, "mmpbsa.in")
    shutil.copy(input_file_src, input_file_dst)
    logger.info(f"Copied GBSA input file from {input_file_src} to {input_file_dst}")
    prep_mmpbsain(logger, input_file_dst, start, stop, dna, phos, ffs)
    logger.debug(f"Prepared GBSA input file: {input_file_dst}")
    logger.debug(f"Computing GBSA with selections: {index_group1} and {index_group2}")


    # Run MMPBSA calculation
    cmd = ["mpirun", "--use-hwthread-cpus", 
           "-np", str(cores),
           "gmx_MMPBSA", 
           "-i", Path(input_file_dst).resolve(), 
           "-cs", Path(tprfile).resolve(), 
           "-ci", Path(index).resolve(), 
           "-cg", index_group1, index_group2,
           "-ct", Path(xtcfile).resolve(), 
           "-cp", Path(topfile).resolve(), 
           "-nogui", 
           "-o", "FINAL_RESULTS_MMPBSA.dat",
           "-do", "FINAL_DECOMP_MMPBSA.dat",
           "-eo", "FINAL_ENERGY.csv", 
           "-deo", "FINAL_ENERGY_DECOMP.csv"
           ]
    
    subprocess.run(cmd, cwd=outdir, check=True, text=True)

if __name__ == "__main__":
    main()