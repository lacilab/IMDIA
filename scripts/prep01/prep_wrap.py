#!/usr/bin/env python
# -*- coding: utf-8 -*-

from datetime import datetime
import subprocess
import sys
from scripts.utils.helpers import gen_flags, setup_logger
from scripts.prep01 import prep

def main(datadir: str, filename: str,
    dna : bool = False, phos : bool = False,
    range_a : str = None, range_b : str = None,
    ffs : str = None):

    args = locals()
    args.pop('datadir', None)
    args.pop('filename', None)
    flags = gen_flags(args)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    LOGFILE = f"{datadir}/logs/step01_prep_{timestamp}.log"
    logger = setup_logger("prep", LOGFILE)
    logger.debug("Running preparation step.")

    # cmd = [sys.executable,
    #     "-m",
    #     "scripts.prep01.prep",
    #     "--logfile", LOGFILE,
    #     "--datadir", datadir,
    #     "--filename", filename,
    # ] + flags

    prep.main(logfile=LOGFILE, datadir=datadir, filename=filename,
        dna=dna, phos=phos, range_a=range_a, range_b=range_b, ffs=ffs)

    # remove empty strings from cmd before joining
    # cmd = [c for c in cmd if c != ""]
    # subprocess.run(cmd, check=True)