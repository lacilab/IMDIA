import cmd
import os
import subprocess
from utils.helpers import gen_flags

def main(gmxdir: str, repn: int, 
    dna = False, phos = False,
    ntmpi = None, ntomp = None, ngpu = None, gpulist = None, 
    em = None, nvt = None, npt = None, md = None, 
    ):

    args = locals()
    args.pop('gmxdir')
    flags = gen_flags(args)

    here = os.path.dirname(os.path.abspath(__file__))
    sh_path = os.path.join(here, "md.py")
    
    cmd = [
        sh_path,
        "--gmxdir", gmxdir
    ] + flags

    print("Running command for md script:")
    # remove empty strings from cmd before joining
    cmd = [c for c in cmd if c != ""]
    print(" ".join(cmd))
    subprocess.call(cmd)