#!/usr/bin/env python3

import sys
import os
import re
import pickle
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import numpy as np
import pandas as pd
from scripts.utils.helpers import setup_logger


def main(logfile: str = None, universe: str = None, 
         tprfile: str = None, xtcfile: str = None, 
         outdir : str = None, hbonds_mda_parsed : str = None,
         start: int = 0, stop: int = None,
         group1: str = None, group2: str = None, dna : bool = False):
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
    logger.debug("Starting hydrogen bond analysis with MDAnalysis.")

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
        logger.info(f"Created output directory: {outdir}")
    else:
        logger.warning(f"Output directory already exists: {outdir}")

    # Check that group1 and group2 are provided
    if group1 is None or group2 is None:
        logger.error("Both group1 and group2 must be specified for interaction analysis.")
        sys.exit(1)

    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[0]

    u_tpr = mda.Universe(tprfile,xtcfile)

    # Define selections
    selections = []
    for group in [group1, group2]:
        if re.match(r"^[A-Za-z]$", group):
            logger.debug(f"MDA - hbond analysis - Group {group} identified as chain.")
            selections.append(u.select_atoms(f"segid {group}"))
        elif group == "DNA" and dna:
            logger.debug(f"MDA - hbond analysis - Group {group} identified as DNA.")
            selections.append(u.select_atoms("nucleic")) 
        elif group == "Protein":
            logger.debug(f"MDA - hbond analysis - Group {group} identified as Protein.")
            protein_sel = u.select_atoms("protein")
            # Get segids of the protein, to avoid leaving out non-standard residues like PTR
            prot_segids = [seg.segid for seg in protein_sel.segments]
            logger.debug(f"Protein segids identified: {prot_segids}")
            selections.append(u.select_atoms("segid " + " ".join(prot_segids)))

    # Convert AtomGroup indices to selection strings
    sele_index = []
    for sel in selections:
        sele_index.append("index " + " ".join(map(str, sel.indices)))

    # Define hbond analysis parameters
    hbonds_group_a_to_group_b = HBA(
        universe=u_tpr,
        between=[sele_index],
        d_h_a_angle_cutoff=120,
        d_a_cutoff=3.5
    )
    # Run analysis, starting at frame 1000 until 2000 (end)
    hbonds_group_a_to_group_b.run(start=start, stop=stop, verbose=True)
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