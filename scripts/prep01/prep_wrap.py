#!/usr/bin/env python
# -*- coding: utf-8 -*-

from datetime import datetime
from pathlib import Path
from scripts.utils.helpers import setup_logger
from scripts.prep01 import prep
from scripts.pipeconfig import Pipeline

def main(pipe: Pipeline):
    """ 
    Run preparation module to prepare a PDB structure for MD simulations.

    Steps:
        1.Amber preparation: setup forcefields, add hydrogens.
        2.Conversion from Amber to GROMACS.
        3.Extract .itp files from .top topology file.
        4.Remove ions added with Amber.
        5.Add ions with GROMACS to neutralize the system and reach target concentration.
    """

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logdir = pipe.paths.mutdir / "logs"
    logdir.mkdir(parents=True, exist_ok=True)
    LOGFILE = logdir / f"step01_prep_{timestamp}.log"
    logger = setup_logger("prep", LOGFILE)
    logger.debug("Running preparation step.")

    prep.main(logfile=LOGFILE, datadir=pipe.paths.mutdir, filename=pipe.paths.pdbfile,
        dna=pipe.system.dna, phos=pipe.system.phos, 
        range_a=pipe.system.range_a, range_b=pipe.system.range_b, ffs=pipe.system.ffs)
    