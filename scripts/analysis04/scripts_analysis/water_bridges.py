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
from scripts.utils.parallel import an_perblock

# Custom forcefield atom names
class WaterBridgeAnalysis_OtherFF(WaterBridgeAnalysis):
    DEFAULT_DONORS = {"OtherFF": tuple(set(["O", "H2", "N", "NE", "NH1", "NH2", "ND2", "SG", "NE2", "ND1", 
                                            "NZ", "OG", "OG1", "OH"]))}
    DEFAULT_ACCEPTORS = {"OtherFF": tuple(set(["O", "H2", "OC1", "OC2", "OD1", "OD2", 
                                                "SG", "OE1", "OE2", "ND1", "NE2", "SD", "OG", "OG1", "OH"]))}

def water_bridge_factory(u, start, stop, sel1, sel2, distance, angle):
    wba = WaterBridgeAnalysis_OtherFF(
            universe=u,
            selection1=sel1,
            selection2=sel2,
            water_selection="resname WAT",  
            forcefield="OtherFF",
            distance=distance,  
            angle=angle,     
    )

    # Run the analysis over the trajectory
    wba.run(start=start, stop=stop, verbose=True)
    return wba.generate_table()


def main(logfile : str = None, universe : str = None, tprfile: str = None, xtcfile: str = None, outcsv : str = None, cores : int = 1,
         group1: str = None, group2: str = None, 
         start: int = None, stop: int = None, dna : bool = False):
    """
    Analyze water-mediated bridges between target groups in an MD trajectory.

    Arguments:
        logfile   : Path to the log file.
        universe  : Path to the Universe object generated with MDAnalysis, containing the processed MD system.
        tprfile   : GROMACS run input file (.tpr) containing topology and run parameters.
        xtcfile   : GROMACS trajectory file (.xtc) containing trajectory data.
        outcsv    : Output .csv file containing water bridge analysis results.

        cores     : Number of CPU cores used for parallel processing.

        group1    : First atom selection group for analysis.
        group2    : Second atom selection group for analysis.

        start     : Starting frame index for analysis.
        stop      : Ending frame index for analysis.

        dna       : True if DNA is present in the system.
    """

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug("Starting water bridge analysis with MDAnalysis.")

    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[0]


    # Define selections
    selections = []
    for group in [group1, group2]:
        if re.match(r"^[A-Za-z]$", group):
            logger.debug(f"MDA - water bridge analysis - Group {group} identified as chain.")
            selections.append(f"segid {group}")
        elif group == "DNA" and dna:
            logger.debug(f"MDA - water bridge analysis - Group {group} identified as DNA.")
            selections.append("nucleic") 
        elif group == "Protein":
            logger.debug(f"MDA - water bridge analysis - Group {group} identified as Protein.")
            protein_sel = u.select_atoms("protein")
            # Get segids of the protein, to avoid leaving out non-standard residues like PTR
            prot_segids = [seg.segid for seg in protein_sel.segments]
            logger.debug(f"Protein segids identified: {prot_segids}")
            selections.append("segid " + " ".join(prot_segids))

    if len(selections) != 2:
        logger.warning("Error in defining selections for water bridge analysis.")

    sel1, sel2 = selections
    distance = 3.5
    angle = 120

    if cores > 1:
        logger.debug(f"MDA - water bridge analysis - Running parallel Dask with {cores} cores")
        table = an_perblock(logger, tprfile, xtcfile, cores, water_bridge_factory, start, stop, sel1, sel2, distance, angle)
    else:
        logger.debug(f"MDA - water bridge analysis - Running single-threaded")
        # Initialize WaterBridgeAnalysis with custom forcefield atom names
        wba = WaterBridgeAnalysis_OtherFF(
                universe=u,
                selection1=sel1,
                selection2=sel2,
                water_selection="resname WAT",  
                forcefield="OtherFF",
                distance=distance,  
                angle=angle,     
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
