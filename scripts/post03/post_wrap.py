#!/usr/bin/env python
# -*- coding: utf-8 -*-

from datetime import datetime
import subprocess
import sys
import os
from scripts.utils.helpers import gen_flags, setup_logger

def main(datadir: str, repn: int,
    dna : bool = False, phos : bool = False):

    args = locals()
    args.pop('datadir', None)
    args.pop('repn', None)
    flags = gen_flags(args)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    LOGFILE = f"{datadir}/logs/step03_rep{repn}_post_{timestamp}.log"
    logger = setup_logger("post", LOGFILE)

    repdir = f"{datadir}/rep{repn}/gromacs"
    logger.debug(f"Repdir: {repdir}")

    cmd = [sys.executable,
        "-m",
        "scripts.post03.post",
        "--logfile", LOGFILE,
        "--repdir", repdir,
    ] + flags

    logger.debug("Running post-processing step.")
    # remove empty strings from cmd before joining
    cmd = [c for c in cmd if c != ""]
    subprocess.run(cmd, check=True)