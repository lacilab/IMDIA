#!/usr/bin/env python3

import subprocess
import MDAnalysis as mda
import warnings
import pickle
from pathlib import Path
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')
import os
import sys
from scripts.utils.helpers import setup_logger

def frame2ps(u, frame):
    frame_dt_ps = u.trajectory.dt
    return int(frame * frame_dt_ps)

def main(logfile: str =None, universe: str = None, gmxdir : str = None, outdir : str =None, indexfile : str = None,
         group1: str = None, group2: str = None, start: str = None, stop: str = None):
    """
    Predict hydrogen bonds between two groups using GROMACS.

    Arguments:
        logfile    : Path to the log file.
        universe   : Path to the Universe object generated with MDAnalysis, containing the processed MD system.
        gmxdir     : Path to the directory with the MD results (e.g. rep1/gromacs).
        outdir     : Output directory for results (e.g. rep1/analysis/hbonds).
        indexfile  : GROMACS index file (.ndx) defining atom groups.

        group1     : First atom group for hydrogen bond analysis.
        group2     : Second atom group for hydrogen bond analysis.

        start      : Start frame index for analysis window.
        stop       : End frame ndex for analysis window.
    """
    
    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
        logger.debug(f"Created output directory: {outdir}")
    else:
        logger.warning(f"Output directory already exists: {outdir}")

    # Check that group1 and group2 are provided
    if group1 is None or group2 is None:
        logger.error("Both group1 and group2 must be specified for interaction analysis.")
        sys.exit(1)

    if start is not None or stop is not None:
        with open(universe, 'rb') as f:
            u = pickle.load(f)
            u.trajectory[0]
    start_option = ["-b", str(frame2ps(u, start))] if start is not None else []
    stop_option = ["-e", str(frame2ps(u, stop))] if stop is not None else []
    if len(start_option) == 0:
        logger.debug("gmx hbond: No start frame specified, using the beginning of the trajectory.")
    else:
        logger.debug(f"gmx hbond: Using start frame {start} corresponding to {start_option} ps")
    if len(stop_option) == 0:
        logger.debug("gmx hbond: No stop frame specified, using the end of the trajectory.")
    else:
        logger.debug(f"gmx hbond: Using stop frame {stop} corresponding to {stop_option} ps")

    outfile = Path(outdir) / "hbonds_gmx_num.xvg"

    # Calculate hydrogen bonds between selections
    cmd = ["gmx", "hbond", 
           "-f", f"{gmxdir}/md_post.xtc", 
           "-s", f"{gmxdir}/md.tpr", 
           "-n", indexfile, 
           *start_option,
           *stop_option,
           "-o", f"{outdir}/hbonds_gmx.ndx", 
           "-pf", "yes", 
           "-num", outfile]
    
    logger.debug(f"Running gmx hydrogen bonds analysis between groups: {group1} and {group2}")
    process = subprocess.run(cmd, input=f"{group1}\n{group2}\n", text=True, check=True)

    if process.returncode != 0:
        logger.error(f"gmx hbond failed (exit code {process.returncode}):\n{process.stderr}")
    else:
        logger.debug(f"gmx hbond output:\n{process.stdout}\n{process.stderr}")

    if outfile.is_file():
        logger.info(f"Hydrogen bonds analysis completed with GROMACS. Results saved in: {outfile}")
    else:
        logger.warning(f"Hydrogen bonds analysis with GROMACS did not produce {outfile}. "
                    f"Check log {logfile}.")    


if __name__ == "__main__":
    main()