#!/usr/bin/env python3

import multiprocessing
import os
import shutil
import subprocess
from pathlib import Path
from scripts.utils.helpers import setup_logger
from scripts.md02.group_index import main as group_index

def main(logfile: str = None, scriptdir: str = None, 
         gmxdir: str = None, outdir: str = None, tprfile: str = None,
         xtcfile: str = None, topfile: str = None,
         dna: bool = False, phos: bool = False):
    # Initialize logger
    logger = setup_logger("analysis", logfile)

    # Create output directory if it doesn't exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        logger.info(f"Created output directory: {outdir}")
    else:
        logger.warning(f"Directory already exists: {outdir}")

    # Copy input file
    input_file_src = os.path.join(scriptdir, "mmpbsa.in")
    input_file_dst = os.path.join(outdir, "mmpbsa.in")
    shutil.copy(input_file_src, input_file_dst)

    # Prepare index file
    protein_group = group_index(gmxdir=gmxdir, logfile=logfile, logname='analysis', structurefile=f"{gmxdir}/md_sys.pdb", dna=dna, phos=phos)

    if dna:
        sele1 = 'DNA'
        sele2 = protein_group
        logger.info(f"Preparing MMPBSA index file for protein-DNA system, with index groups: {sele1}, {sele2}")
    else:
        sele1 = 'Chain_A'
        sele2 = 'Chain_B'
        logger.info(f"Preparing MMPBSA index file for protein-protein system, with index groups: {sele1}, {sele2}")

    tprfile_abs = Path(tprfile).resolve()
    xtcfile_abs = Path(xtcfile).resolve()
    topfile_abs = Path(topfile).resolve()
    gmxdir_abs = Path(gmxdir).resolve()

    # Auto-detect NTMPI and NTOMP if not provided
    SLURM_NTASKS = os.environ.get('SLURM_NTASKS')
    SLURM_CPUS_PER_TASK = os.environ.get('SLURM_CPUS_PER_TASK')
    if SLURM_NTASKS and SLURM_CPUS_PER_TASK:
        NP = int(SLURM_NTASKS)
    else:
        # fallback for local machine
        NP = multiprocessing.cpu_count()

    # Run MMPBSA calculation
    cmd = ["mpirun", "--use-hwthread-cpus", 
           "-np", NP, "gmx_MMPBSA", 
           "-i", "mmpbsa.in", 
           "-cs", tprfile_abs, 
           "-ci", f"{gmxdir_abs}/index.ndx", 
           "-cg", sele1, sele2,
           "-ct", xtcfile_abs, 
           "-cp", topfile_abs, 
           "-nogui", 
           "-o", "FINAL_RESULTS_MMPBSA.dat",
           "-do", "FINAL_DECOMP_MMPBSA.dat",
           "-eo", "FINAL_ENERGY.csv", 
           "-deo", "FINAL_ENERGY_DECOMP.csv"
           ]
    
    subprocess.run(cmd, cwd=outdir, check=True, text=True)

if __name__ == "__main__":
    main()