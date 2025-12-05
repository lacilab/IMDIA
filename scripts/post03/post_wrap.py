#!/usr/bin/env python
# -*- coding: utf-8 -*-

from datetime import datetime
import subprocess
import sys
from scripts.utils.helpers import gen_flags, setup_logger

def main(datadir: str, gmxdir: str,
    dna = False, phos = False):

    args = locals()
    args.pop('datadir', None)
    args.pop('gmxdir', None)
    flags = gen_flags(args)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    LOGFILE = f"{datadir}/logs/step03_post_{timestamp}.log"
    logger = setup_logger("post", LOGFILE)

    cmd = [sys.executable,
        "-m",
        "scripts.post03.post",
        "--logfile", LOGFILE,
        "--gmxdir", gmxdir,
    ] + flags

    logger.debug("Running post-processing step.")
    # remove empty strings from cmd before joining
    cmd = [c for c in cmd if c != ""]
    subprocess.run(cmd, check=True)