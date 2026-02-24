#!/usr/bin/env python3

from os import path
import os
from pathlib import Path
import shutil
import subprocess
from scripts.utils.helpers import setup_logger

def main(logfile : str = None, tprfile : Path = None, xtcfile : Path = None, scriptdir: Path = None, outdir : Path = None):
    ##### Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug(f"The working directory is: {os.getcwd()}")
    logger.debug(f"Input tprfile: {tprfile}")
    logger.debug(f"Input xtcfile: {xtcfile}")
    logger.debug(f"Input scriptdir: {scriptdir}")
    logger.debug(f"Input outdir: {outdir}")

    # Check if directories and extensions
    if not path.isfile(tprfile):
        logger.error(f"tpr file not found: {tprfile}")
        return
    elif tprfile.suffix != ".tpr":
        logger.error(f"tpr file does not have .tpr extension: {tprfile}")
        return
    if not path.isfile(xtcfile):
        logger.error(f"xtc file not found: {xtcfile}")
        return
    elif xtcfile.suffix != ".xtc":
        logger.error(f"xtc file does not have .xtc extension: {xtcfile}")
        return
    if not path.isdir(scriptdir):
        logger.error(f"Script directory not found: {scriptdir}")
        return
    if not path.isdir(outdir):
        logger.info(f"Output directory not found, creating: {outdir}")
        os.makedirs(outdir, exist_ok=True)

    shutil.copyfile(f"{scriptdir}/porcupine.pml", f"{outdir}/porcupine.pml")
    if path.isfile(f"{outdir}/porcupine.pml"):
        logger.info(f"Copied porcupine.pml to: {outdir}/porcupine.pml")
    shutil.copyfile(f"{scriptdir}/modevectors.py", f"{outdir}/modevectors.py")
    if path.isfile(f"{outdir}/modevectors.py"):
        logger.info(f"Copied modevectors.py to: {outdir}/modevectors.py")

    # Generate extreme1.gro using gmx anaeig
    cmd = ["gmx", "anaeig", 
           "-v", "eigenvec.trr", 
           "-s", tprfile.resolve(), 
           "-f", xtcfile.resolve(), 
           "-b", "0", 
           "-first", "1", "-last", "1", 
           "-extr", "extreme1.gro"
    ]
    subprocess.run(cmd, cwd=outdir, check=True, text=True, input="C-alpha\nC-alpha\n")
    logger.info(f"Executed gmx anaeig in directory: {outdir}")
    if path.isfile(f"{outdir}/extreme1.gro"):
        logger.info(f"Generated extreme1.gro at: {outdir}/extreme1.gro")
    
    # Visualise using modevectors.py, no GUI
    cmd = [
        "pymol", "-cq", "porcupine.pml"
    ]
    subprocess.run(cmd, cwd=outdir, check=True, text=True)
    if path.isfile(f"{outdir}/porcupine.png"):
        logger.info(f"Porcupine plot generated at: {outdir}/porcupine.png")

    # Clean temporary files
    os.remove(f"{outdir}/porcupine.pml")
    os.remove(f"{outdir}/modevectors.py")

if __name__ == "__main__":
    main()



