#!/usr/bin/env python
# -*- coding: utf-8 -*-

from datetime import datetime
from scripts.utils.helpers import setup_logger
from scripts.pipeconfig import Pipeline
from scripts.analysis04 import analysis

def main(pipe: Pipeline):
    """ Run analysis module. """
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logdir = pipe.paths.mutdir / "logs"
    logdir.mkdir(parents=True, exist_ok=True)  # Ensure dir exists
    LOGFILE = logdir / f"step04_analysis_{timestamp}.log"
    logger = setup_logger("analysis", LOGFILE)
    logger.debug(f"Running analysis step with rules: {pipe.rules} and unlock: {pipe.unlock}.")

    analysis.main(logfile=LOGFILE, datadir=pipe.paths.datadir, configfile=pipe.paths.paramspath, pdbpath=pipe.pdbpath,
                  postpdb=pipe.paths.postpdb, cores=pipe.run.snakecores, rules=pipe.rules, unlock=pipe.unlock,
                  flags={"protein": pipe.system.protein, "dna": pipe.system.dna, "rna": pipe.system.rna, "phos": pipe.system.phos})

