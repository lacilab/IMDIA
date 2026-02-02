#!/usr/bin/env python3

import os
import sys
import re
import pickle
import argparse
from venv import logger
import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import WaterBridgeAnalysis
from scripts.utils.helpers import setup_logger


def main(logfile : str = None, universe : str = None, outcsv : str = None, 
         group1: str = None, group2: str = None, 
         start: int = None, stop: int = None, dna : bool = False):
    # Arguments
    # parser = argparse.ArgumentParser(description="Analyze water bridges between protein and nucleic acids.")
    # parser.add_argument("-u", "--universe", help="Path to the Universe object")
    # parser.add_argument("--dna", action='store_true', dest='with_dna', help="DNA is present in the system")
    # parser.add_argument("-o", "--wb_csv", help="Output CSV file with water bridges data.")

    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit(0)
 
    # args = parser.parse_args()
    # universe = args.universe
    # wb_csv = args.wb_csv

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug("Starting water bridge analysis with MDAnalysis.")

    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[0]

    # Custom forcefield atom names
    class WaterBridgeAnalysis_OtherFF(WaterBridgeAnalysis):
        DEFAULT_DONORS = {"OtherFF": tuple(set(["O", "H2", "N", "NE", "NH1", "NH2", "ND2", "SG", "NE2", "ND1", 
                                                "NZ", "OG", "OG1", "OH"]))}
        DEFAULT_ACCEPTORS = {"OtherFF": tuple(set(["O", "H2", "OC1", "OC2", "OD1", "OD2", 
                                                    "SG", "OE1", "OE2", "ND1", "NE2", "SD", "OG", "OG1", "OH"]))}

    # Define selections
    selections = []
    for group in [group1, group2]:
        if re.match(r"^[A-Za-z]$", group):
            logger.debug(f"MDA - contact analysis - Group {group} identified as chain.")
            selections.append(f"segid {group}")
        elif group == "DNA" and dna:
            logger.debug(f"MDA - contact analysis - Group {group} identified as DNA.")
            selections.append("nucleic") 
        elif group == "Protein":
            logger.debug(f"MDA - contact analysis - Group {group} identified as Protein.")
            protein_sel = u.select_atoms("protein")
            # Get segids of the protein, to avoid leaving out non-standard residues like PTR
            prot_segids = [seg.segid for seg in protein_sel.segments]
            logger.debug(f"Protein segids identified: {prot_segids}")
            selections.append("segid " + " ".join(prot_segids))

    if len(selections) != 2:
        logger.warning("Error in defining selections for contact analysis.")

    sel1, sel2 = selections

    # Initialize WaterBridgeAnalysis with custom forcefield atom names
    wba = WaterBridgeAnalysis_OtherFF(
            universe=u,
            selection1=sel1,
            selection2=sel2,
            water_selection="resname WAT",  
            forcefield="OtherFF",
            distance=3.5,  
            angle=120,     
    )


    # Run the analysis over the trajectory
    wba.run(start=start, stop=stop, verbose=True)

    # Generate the table
    table = wba.generate_table()

    # Convert and save
    df = pd.DataFrame(table)
    df.to_csv(outcsv, index=False)

    # Check output file
    if os.path.exists(outcsv):
        logger.info(f"Water bridges data saved to: {outcsv}")
    else:
        logger.warning(f"Failed to save water bridges data to: {outcsv}")
