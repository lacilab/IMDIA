import os
import subprocess
import sys
import datetime
from scripts.utils.helpers import gen_flags

def main(datadir: str, gmxdir: str, repn: int, 
    dna = False, phos = False,
    ntmpi = None, ntomp = None, ngpu = None, gpulist = None, 
    em = None, nvt = None, npt = None, md = None, 
    ):

    args = locals()
    args.pop('datadir')
    args.pop('gmxdir')
    flags = gen_flags(args)
    
    cmd = [
        sys.executable,
        "-m",
        "scripts.md02.md",
        "--gmxdir", gmxdir
    ] + flags

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    LOGFILE = os.path.join(datadir, "logs/", f"rep{repn}_step02_md_{timestamp}.log")

    print("Running command for md script:")
    # remove empty strings from cmd before joining
    cmd = [x for x in cmd if x.strip()]
    with open(LOGFILE, "w") as log:
        subprocess.run(cmd, stdout=log, stderr=log, check=True)

