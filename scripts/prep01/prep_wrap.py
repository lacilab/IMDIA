import os
import subprocess
from utils.helpers import gen_flags 

def main(datadir: str, filename: str,
    dna = False, phos = False,
    range_a = None, range_b = None,
    ffs = None):

    args = locals()
    args.pop(['datadir', 'filename'])
    flags = gen_flags(args)

    here = os.path.dirname(os.path.abspath(__file__))
    sh_path = os.path.join(here, "prep_00_all.sh")

    cmd = [
        sh_path,
        datadir,
        filename,
    ] + flags

    print("Running command for prep script:")
    # remove empty strings from cmd before joining
    cmd = [c for c in cmd if c != ""]
    print(" ".join(cmd))
    subprocess.call(cmd)