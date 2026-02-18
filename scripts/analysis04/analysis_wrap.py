#!/usr/bin/env python
# -*- coding: utf-8 -*-

from datetime import datetime
from scripts.utils.helpers import setup_logger
from scripts.pipeconfig import Pipeline
from scripts.analysis04 import analysis

def main(pipe: Pipeline, rules: list[str] = None, unlock: bool = False):

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logdir = pipe.paths.mutdir / "logs"
    logdir.mkdir(parents=True, exist_ok=True)  # Ensure dir exists
    LOGFILE = logdir / f"step04_analysis_{timestamp}.log"
    logger = setup_logger("analysis", LOGFILE)
    logger.debug(f"Running analysis step with rules: {rules} and unlock: {unlock}.")

    analysis.main(logfile=LOGFILE, datadir=pipe.paths.datadir, configfile=pipe.paths.configfile, 
                  cores=pipe.run.cores, rules=rules, unlock=unlock)

