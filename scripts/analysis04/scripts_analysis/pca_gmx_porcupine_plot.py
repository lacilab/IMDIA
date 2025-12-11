#!/usr/bin/env python3

from os import path
import os
from pathlib import Path
import shutil
import subprocess
from scripts.utils.helpers import setup_logger

def main(logfile : str = None, tprfile : str = None, xtcfile : str = None, scriptdir: str = None, outdir : str = None):
    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    # TODO check if directories and files exist
    # TODO check that the files have the correct extensions

    shutil.copyfile(f"{scriptdir}/porcupine.pml", f"{outdir}/porcupine.pml")
    if path.isfile(f"{outdir}/porcupine.pml"):
        logger.info(f"Copied porcupine.pml to: {outdir}/porcupine.pml")
    shutil.copyfile(f"{scriptdir}/modevectors.py", f"{outdir}/modevectors.py")
    if path.isfile(f"{outdir}/modevectors.py"):
        logger.info(f"Copied modevectors.py to: {outdir}/modevectors.py")

    tprfile_abs = Path(tprfile).resolve()
    xtcfile_abs = Path(xtcfile).resolve()

    # Generate extreme1.gro using gmx anaeig
    cmd = ["gmx", "anaeig", 
           "-v", "eigenvec.trr", 
           "-s", tprfile_abs, 
           "-f", xtcfile_abs, 
           "-b", "0", "-e", "200001", 
           "-first", "1", "-last", "1", 
           "-extr", "extreme1.gro", "-nframes", "2001"
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



