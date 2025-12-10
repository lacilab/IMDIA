#!/usr/bin/env python3

import subprocess
import MDAnalysis as mda
import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')
import os
from scripts.utils.helpers import setup_logger
from scripts.md02.group_index import main as group_index

def main(logfile: str =None, gmxdir : str = None, outdir : str =None, dna : bool = False, phos : bool = False):
    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
        logger.info(f"Created output directory: {outdir}")
    else:
        logger.warning(f"Output directory already exists: {outdir}")
    
    # Prepare index file
    protein_group = group_index(gmxdir=gmxdir, logfile=logfile, logname='analysis', structurefile=f"{gmxdir}/md_sys.pdb", dna=dna, phos=phos)

    # Calculate hydrogen bonds between selections
    cmd = ["gmx", "hbond", 
           "-f", f"{gmxdir}/md_noPBC_mol_com_cluster_fit.xtc", 
           "-s", f"{gmxdir}/md.tpr", 
           "-n", f"{gmxdir}/index.ndx", 
           "-b", "100000",
           "-o", f"{outdir}/hbonds_gmx.ndx", 
           "-pf", "yes", 
           "-num", f"{outdir}/hbonds_gmx_num.xvg"]
    
    if dna:
        input_sele = f"""
        {protein_group}
        DNA
        """
    else:
        input_sele = """
        Chain_A
        Chain_B
        """

    subprocess.run(cmd, input=input_sele, text=True, check=True)
    logger.info(f"Hydrogen bonds analysis completed with GROMACS. Results saved in: {outdir}")
    


if __name__ == "__main__":
    main()