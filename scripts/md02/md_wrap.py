#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
from datetime import datetime
from scripts.utils.helpers import setup_logger
from scripts.utils.group_index import group_index
from scripts.pipeconfig import Pipeline
from scripts.md02 import md


def main(pipe: Pipeline, steps: list):
    """ Perform MD simulation steps: minimization, NVT, NPT, and production MD. """ 
    
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
        if "md" in steps:
            logger.debug(f"Running full MD for replicate {n}.")
        else:
            logger.debug(f"Running specified MD steps {steps} for replicate {n}.")
        flag_nb_gpu, flag_pme_gpu, flag_gpu_tasks = md.setup(logfile=LOGFILE, 
                                                             gmxdir=pipe.paths.gmxdir, outdir=outdir, repn=n,
                                                             ntmpi=pipe.run.ntmpi, ntomp=pipe.run.ntomp, 
                                                             nb=pipe.run.nb, pme=pipe.run.pme, npme=pipe.run.npme,
                                                             ngpu=pipe.run.ngpu, gputasks=pipe.run.gputasks,
                                                             em=pipe.run.em, nvt=pipe.run.nvt, npt=pipe.run.npt, md=pipe.run.md)
        if "md_em" in steps or "md" in steps:
            md.minimization(logfile=LOGFILE, gmxdir=pipe.paths.gmxdir, outdir=outdir, flag_nb_gpu=flag_nb_gpu)
        if "md_eq" in steps or "md" in steps:
            group_index(gmxdir=outdir, logfile=LOGFILE, logname='md', structurefile=outdir / "em.gro", dna=pipe.system.dna, phos=pipe.system.phos)
            md.equilibration(logfile=LOGFILE, gmxdir=pipe.paths.gmxdir, outdir=outdir, flag_nb_gpu=flag_nb_gpu)
        if "md_prod" in steps or "md" in steps:
            print(f"LOGFILE type: {type(LOGFILE)}, value: {LOGFILE}")
            print(f"gmxdir type: {type(pipe.paths.gmxdir)}, value: {pipe.paths.gmxdir}")
            print(f"outdir type: {type(outdir)}, value: {outdir}")
            print(f"md_time type: {type(pipe.run.md_time)}, value: {pipe.run.md_time}")
            group_index(gmxdir=outdir, logfile=LOGFILE, logname='md', structurefile=outdir / "em.gro", dna=pipe.system.dna, phos=pipe.system.phos)
            md.production(logfile=LOGFILE, gmxdir=pipe.paths.gmxdir, outdir=outdir, md_time=str(pipe.run.md_time), ntmpi=pipe.run.ntmpi, ntomp=pipe.run.ntomp, flag_nb_gpu=flag_nb_gpu, flag_pme_gpu=flag_pme_gpu, flag_gpu_tasks=flag_gpu_tasks)