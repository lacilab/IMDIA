#!/usr/bin/env python
# -*- coding: utf-8 -*-

from datetime import datetime
from pathlib import Path
from scripts.utils.helpers import setup_logger
from scripts.diag00 import diag
from scripts.pipeconfig import Pipeline

def main(pipe: Pipeline):
    """
    Run diagnosis module to check a PDB structure before preparation.

    Steps:
        1.Gap check: chain breaks, backbone distortions, missing atoms, numbering jumps.
    """

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logdir = pipe.paths.mutdir / "logs"
    logdir.mkdir(parents=True, exist_ok=True)
    LOGFILE = logdir / f"step00_diag_{timestamp}.log"
    logger = setup_logger("diag", LOGFILE)
    logger.info("Running diagnosis step.")

    diag.main(logfile=LOGFILE, datadir=pipe.paths.mutdir, filename=pipe.paths.pdbfile,
        stop_on_gaps=pipe.system.stop_on_gaps)
