#!/usr/bin/env python3

import os
import shutil
import subprocess
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

    # Run MMPBSA calculation
    cmd = ["mpirun", "--use-hwthread-cpus", 
           "-np", "4", "gmx_MMPBSA", 
           "-i", f"{outdir}/mmpbsa.in", 
           "-cs", tprfile, 
           "-ci", f"{gmxdir}/index.ndx", 
           "-cg", sele1, sele2,
           "-ct", xtcfile, 
           "-cp", topfile, 
           "-nogui", 
           "-prefix", f"{outdir}/_GMXMMPBSA_",
           "-o", f"{outdir}/FINAL_RESULTS_MMPBSA.dat",
           "-do", f"{outdir}/FINAL_DECOMP_MMPBSA.dat",
           "-eo", f"{outdir}/FINAL_ENERGY.csv", 
           "-deo", f"{outdir}/FINAL_ENERGY_DECOMP.csv"
           ]
    subprocess.run(cmd, check=True, text=True)

if __name__ == "__main__":
    main()