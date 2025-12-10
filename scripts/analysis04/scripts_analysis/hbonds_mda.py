#!/usr/bin/env python3

import sys
import os
import pickle
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import numpy as np
import pandas as pd
from scripts.utils.helpers import setup_logger


def main(logfile: str = None, universe: str = None, 
         tprfile: str = None, xtcfile: str = None, 
         outdir : str = None, hbonds_mda_parsed : str = None,
         dna : bool = False):
    # # Arguments
    # parser = argparse.ArgumentParser(description="Parse MDAnalysis hydrogen bonds output into a CSV file.")
    # parser.add_argument("-u", "--universe", help="MDAnalysis universe file.")
    # parser.add_argument("-t", "--tpr", help="Path to the TPR file.")
    # parser.add_argument("-x", "--xtc", help="Path to the XTC file.")
    # parser.add_argument("--dna", action='store_true', dest='with_dna', help="DNA is present in the system")
    # parser.add_argument("-d", "--outdir", help="Output directory.")
    # parser.add_argument("-o", "--outcsv", help="Output path to save the parsed CSV file.")
# 
    # args = parser.parse_args()
    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit(0)
# 
    # universe = args.universe
    # tprfile = args.tpr
    # xtcfile = args.xtc
    # outdir = args.outdir
    # hbonds_mda_parsed = args.outcsv

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    os.makedirs(outdir, exist_ok=True)

    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[0]

    u_tpr = mda.Universe(tprfile,xtcfile)

    # Define groups
    if dna:
        group_a = u.select_atoms("segid A B") 
        group_b = u.select_atoms("segid C D")
        logger.debug(f"Protein atoms, chains A and B: {len(group_a)}")
        logger.debug(f"DNA atom, chains C and D: {len(group_b)}")
    else:
        group_a = u.select_atoms("segid A") 
        group_b = u.select_atoms("segid B")
        logger.debug(f"Protein chain A atoms: {len(group_a)}")
        logger.debug(f"Protein chain B atoms: {len(group_b)}")

    # Convert AtomGroup indices to selection strings
    group_a_indices = " ".join(map(str, group_a.indices))
    group_b_indices = " ".join(map(str, group_b.indices))

    # Create selection strings
    group_a_sel = f"index {group_a_indices}"
    group_b_sel = f"index {group_b_indices}"

    # Define hbond analysis parameters
    hbonds_group_a_to_group_b = HBA(
        universe=u_tpr,
        between=[
            [group_a_sel, group_b_sel], # for protein-water hbonds
        ],
        d_h_a_angle_cutoff=120,
        d_a_cutoff=3.5
    )
    # Run analysis, starting at frame 1000 until 2000 (end)
    hbonds_group_a_to_group_b.run(start=1000, stop=2001, verbose=True)
    # Extract results
    hbonds_array = hbonds_group_a_to_group_b.results.hbonds

    # Name columns and save raw data
    df = pd.DataFrame(hbonds_array, columns=["frame", "donor", "hydrogen", "acceptor", "distance", "angle"])
    # Adjust to 1-index
    df['donor'] = (df['donor'] + 1).astype('int')
    df['hydrogen'] = (df['hydrogen'] + 1).astype('int')
    df['acceptor'] = (df['acceptor'] + 1).astype('int')

    # Save
    df.to_csv(hbonds_mda_parsed, index=False)
    logger.info(f"Parsed hydrogen bonds saved to {hbonds_mda_parsed}")


if __name__ == "__main__":
    main()