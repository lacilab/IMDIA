#!/usr/bin/env python
# -*- coding: utf-8 -*-

from datetime import datetime
import subprocess
import sys
from scripts.utils.helpers import setup_logger

def main(datadir: str, repn: int, configfile: str, cores: int = None, rules: list = None):

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    LOGFILE = f"{datadir}/logs/step04_analysis_{timestamp}.log"
    logger = setup_logger("analysis", LOGFILE)

    flag_cores = "--cores" if cores is not None else ""
    value_cores = str(cores) if cores is not None else ""

    flag_rules = "--rules" if rules is not None else ""
    value_rules = " ".join(rules) if rules is not None else ""

    cmd = [sys.executable,
        "-m",
        "scripts.analysis04.analysis",
        "--logfile", LOGFILE,
        "--configfile", configfile,
        flag_cores, value_cores,
        flag_rules, value_rules
    ]

    logger.debug("Running analysis step.")
    # remove empty strings from cmd before joining
    cmd = [c for c in cmd if c != ""]
    subprocess.run(cmd, check=True)

