#!/usr/bin/env python3

import sys
import pickle
import argparse
import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import WaterBridgeAnalysis
from scripts.utils.helpers import setup_logger


def main(logfile : str = None, universe : str = None, outcsv : str = None, with_dna : bool = False):
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

    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[0]


    # Custom forcefield atom names
    class WaterBridgeAnalysis_OtherFF(WaterBridgeAnalysis):
        DEFAULT_DONORS = {"OtherFF": tuple(set(["O", "H2", "N", "NE", "NH1", "NH2", "ND2", "SG", "NE2", "ND1", 
                                                "NZ", "OG", "OG1", "OH"]))}
        DEFAULT_ACCEPTORS = {"OtherFF": tuple(set(["O", "H2", "OC1", "OC2", "OD1", "OD2", 
                                                    "SG", "OE1", "OE2", "ND1", "NE2", "SD", "OG", "OG1", "OH"]))}


    if with_dna:
        wba = WaterBridgeAnalysis_OtherFF(
            universe=u,
            selection1='protein or resname PTR',  
            selection2='nucleic',  
            water_selection="resname WAT",  
            forcefield="OtherFF",
            distance=3.5,  
            angle=120,     
        )
    else:
        wba = WaterBridgeAnalysis_OtherFF(
            universe=u,
            selection1='segid A',  
            selection2='segid B',  
            water_selection="resname WAT",  
            forcefield="OtherFF",
            distance=3.5,  
            angle=120,     
        )

    # Run the analysis over the trajectory
    wba.run(start=1000, stop=2001, verbose=True)

    # Generate the table
    table = wba.generate_table()

    # Convert to Pandas DataFrame for easy CSV export
    df = pd.DataFrame(table)

    # Save to CSV
    df.to_csv(outcsv, index=False)
    logger.info(f"Water bridges data saved to: {outcsv}")
