import os
import subprocess
import sys
from scripts.utils.helpers import gen_flags 

def main(datadir: str, filename: str,
    dna = False, phos = False,
    range_a = None, range_b = None,
    ffs = None):

    args = locals()
    args.pop('datadir', None)
    args.pop('filename', None)
    flags = gen_flags(args)

    cmd = [sys.executable,
        "-m",
        "scripts.prep01.prep",
        "--datadir", datadir,
        "--filename", filename,
    ] + flags

    print("Running command for prep script:")
    # remove empty strings from cmd before joining
    cmd = [c for c in cmd if c != ""]
    print(" ".join(cmd))
    subprocess.call(cmd)