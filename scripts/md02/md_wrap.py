#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
from datetime import datetime
from pathlib import Path
from scripts.utils.helpers import gen_flags, setup_logger
from scripts.pipeconfig import Pipeline
from scripts.md02 import md


def main(pipe: Pipeline):

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logdir = pipe.paths.mutdir / "logs"
    logdir.mkdir(parents=True, exist_ok=True)  # Ensure dir exists
    LOGFILE = logdir / f"step02_rep{pipe.system.repn}_md_{timestamp}.log"
    logger = setup_logger("md", LOGFILE)
    logger.debug(f"Running MD step.")

    if not pipe.paths.gmxdir.is_dir():
        logger.error(f"GROMACS directory {pipe.paths.gmxdir} does not exist.")
        sys.exit(1)

    parent_dir = pipe.paths.gmxdir.parent

    if pipe.system.repn == None:
        # List all subdirectories in parent matching rep*
        rep_dirs = [
            d for d in os.listdir(parent_dir)
            if os.path.isdir(os.path.join(parent_dir, d)) and re.match(r'rep\d+', d)
        ]
        # Extract numbers
        rep_numbers = [int(re.search(r'rep(\d+)', d).group(1)) for d in rep_dirs]
        # Sort numerically
        rep_numbers.sort()
        # Determine next replication number
        if not rep_numbers:
            reps = [1]
        else:
            reps = [rep_numbers[-1] + 1]
        logger.debug("Next rep numbers:", reps)
    else:
        reps = pipe.system.repn

    for n in reps:
        outdir = parent_dir / f"rep{n}" / "gromacs"
        logger.debug(f"Starting MD for replicate {n}.")
        md.main(logfile=LOGFILE, gmxdir=pipe.paths.gmxdir, outdir=outdir, repn=n,
            dna=pipe.system.dna, phos=pipe.system.phos,
            ntmpi=pipe.run.ntmpi, ntomp=pipe.run.ntomp, 
            nb=pipe.run.nb, pme=pipe.run.pme, npme=pipe.run.npme,
            ngpu=pipe.run.ngpu, gputasks=pipe.run.gputasks,
            em=pipe.run.em, nvt=pipe.run.nvt, npt=pipe.run.npt, md=pipe.run.md)

