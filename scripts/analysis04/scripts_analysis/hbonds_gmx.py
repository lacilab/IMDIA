#!/usr/bin/env python3

import subprocess
import MDAnalysis as mda
import warnings
import pickle
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')
import os
import sys
from scripts.utils.helpers import setup_logger

def frame2ps(u, frame):
    frame_dt_ps = u.trajectory.dt
    return int(frame * frame_dt_ps)

def main(logfile: str =None, universe: str = None, gmxdir : str = None, outdir : str =None, indexfile : str = None,
         group1: str = None, group2: str = None, start: str = None, stop: str = None, dna : bool = False, phos : bool = False):
    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
        logger.info(f"Created output directory: {outdir}")
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

    # Calculate hydrogen bonds between selections
    cmd = ["gmx", "hbond", 
           "-f", f"{gmxdir}/md_post.xtc", 
           "-s", f"{gmxdir}/md.tpr", 
           "-n", indexfile, 
           *start_option,
           *stop_option,
           "-o", f"{outdir}/hbonds_gmx.ndx", 
           "-pf", "yes", 
           "-num", f"{outdir}/hbonds_gmx_num.xvg"]
    
    logger.debug(f"Running gmx hydrogen bonds analysis between groups: {group1} and {group2}")
    subprocess.run(cmd, input=f"{group1}\n{group2}\n", text=True, check=True)
    logger.info(f"Hydrogen bonds analysis completed with GROMACS. Results saved in: {outdir}")
    


if __name__ == "__main__":
    main()