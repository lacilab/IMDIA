#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import sys
from datetime import datetime
from scripts.utils.helpers import gen_flags, setup_logger


def main(datadir: str, gmxdir: str, repn: int, 
    dna: bool = False, phos: bool = False,
    ntmpi: int = None, ntomp: int = None, nb : str = None, pme :str = None, npme : int = None,
    ngpu : int = None, gpulist : str = None, 
    em : str = None, nvt : str = None, npt : str = None, md : str = None, 
    ):

    args = locals()
    args.pop('datadir')
    args.pop('gmxdir')
    flags = gen_flags(args)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    LOGFILE = f"{datadir}/logs/step02_rep{repn}_md_{timestamp}.log"
    logger = setup_logger("md", LOGFILE)

    cmd = [
        sys.executable,
        "-m",
        "scripts.md02.md",
        "--logfile", LOGFILE,
        "--gmxdir", gmxdir
    ] + flags

    logger.debug("Running command for MD step.")
    # remove empty strings from cmd before joining
    cmd = [x for x in cmd if x.strip()]
    subprocess.run(cmd, check=True)

