#!/usr/bin/env python3

import subprocess
import MDAnalysis as mda
import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')
import os
import sys
from scripts.utils.helpers import setup_logger
from scripts.utils.group_index import group_index_interactions

def main(logfile: str =None, gmxdir : str = None, outdir : str =None, indexfile : str = None,
         group1: str = None, group2: str = None, dna : bool = False, phos : bool = False):
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

    # Calculate hydrogen bonds between selections
    cmd = ["gmx", "hbond", 
           "-f", f"{gmxdir}/md_post.xtc", 
           "-s", f"{gmxdir}/md.tpr", 
           "-n", indexfile, 
           "-b", "100000",
           "-o", f"{outdir}/hbonds_gmx.ndx", 
           "-pf", "yes", 
           "-num", f"{outdir}/hbonds_gmx_num.xvg"]

    input_sele = f"""
        {group1}
        {group2}
        """
    logger.debug(f"Running gmx hydrogen bonds analysis between groups: {group1} and {group2}")
    subprocess.run(cmd, input=input_sele, text=True, check=True)
    logger.info(f"Hydrogen bonds analysis completed with GROMACS. Results saved in: {outdir}")
    


if __name__ == "__main__":
    main()