#!/usr/bin/env python
# -*- coding: utf-8 -*-

from datetime import datetime
import re
import subprocess
import sys
import os
from scripts.utils.helpers import gen_flags, setup_logger
from scripts.pipeconfig import Pipeline
from scripts.post03 import post

def main(pipe: Pipeline):
    """ Run postprocessing module to make MD results ready for analysis and visualization. """

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logdir = pipe.paths.mutdir / "logs"
    logdir.mkdir(parents=True, exist_ok=True)  # Ensure dir exists
    LOGFILE = logdir / f"step03_rep{pipe.system.repn[0]}_post_{timestamp}.log"
    logger = setup_logger("post", LOGFILE)

    logger.debug("Running post-processing step.")

    if pipe.system.repn is None:
        logger.warning("No replication numbers specified in configuration.")
        parent_dir = pipe.paths.gmxdir.parent
        rep_dirs = [
            d for d in os.listdir(parent_dir)
            if os.path.isdir(os.path.join(parent_dir, d)) and re.match(r'rep\d+', d)
        ]
        reps = [int(re.search(r'rep(\d+)', d).group(1)) for d in rep_dirs]
        logger.debug(f"Found replication numbers: {reps}")
    else:
        logger.debug(f"Replication numbers for post-processing: {pipe.system.repn}")
        reps = pipe.system.repn
    
    for n in reps:
        repdir = pipe.paths.mutdir / f"rep{n}" / "gromacs"
        logger.debug(f"Repdir for post-processing: {repdir}")
        post.main(logfile=LOGFILE, repdir=repdir, 
                  dna=pipe.system.dna, phos=pipe.system.phos)
