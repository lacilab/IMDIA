from datetime import datetime
import os
import subprocess
import sys
from scripts.utils.helpers import gen_flags, setup_logger

def main(datadir: str, filename: str,
    dna = False, phos = False,
    range_a = None, range_b = None,
    ffs = None):

    args = locals()
    args.pop('datadir', None)
    args.pop('filename', None)
    flags = gen_flags(args)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    LOGFILE = f"{datadir}/logs/step01_prep_{timestamp}.log"
    logger = setup_logger("prep", LOGFILE)

    cmd = [sys.executable,
        "-m",
        "scripts.prep01.prep",
        "--logfile", LOGFILE,
        "--datadir", datadir,
        "--filename", filename,
    ] + flags

    logger.debug("Running preparation step.")
    # remove empty strings from cmd before joining
    cmd = [c for c in cmd if c != ""]
    subprocess.run(cmd, check=True)