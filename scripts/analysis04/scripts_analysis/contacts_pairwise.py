#!/usr/bin/env python3

import sys
import re
import pickle
from MDAnalysis.lib.distances import distance_array
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scripts.utils.helpers import setup_logger
from scripts.utils.parallel import an_perblock


##### HELPER FUNCTIONS #####

##### Distance by residue pairs
def contact_pairs(u, start_frame, stop_frame, cutoff, sel1, sel2):
    # List to store contact data
    contact_data = []

    # Loop over trajectory
    for ts in u.trajectory[start_frame:stop_frame]:
        frame = ts.frame
        # Calculate distance matrix between sel1 and sel2
        dist_matrix = distance_array(sel1.positions, sel2.positions, box=u.dimensions)
        
        # Find pairs within cutoff
        close_pairs = np.where(dist_matrix <= cutoff)
        
        # Collect contact details
        for i, j in zip(close_pairs[0], close_pairs[1]):
            atom1 = sel1[i]
            atom2 = sel2[j]
            contact = [frame,
                f"{atom1.resid}", 
                f"{atom1.resname}", 
                f"{atom1.segid}",
                f"{atom2.resid}", 
                f"{atom2.resname}",
                f"{atom2.segid}",
                dist_matrix[i, j]
            ]
            contact_data.append(contact)
    return contact_data


##### MAIN FUNCTION #####

def main(logfile : str = None, universe: str = None, 
         tprfile: str = None, xtcfile: str = None,
         cores : int = 1,
         group1 : str = None, group2 : str = None,
         start: int = None, stop: int = None,
         dna : bool = False, outfiles : list = []):
    # Arguments
    # parser = argparse.ArgumentParser(description="Perform contact analysis with MDAnalysis.")
    # parser.add_argument("-u", "--universe", help="Path to the Universe object")
    # parser.add_argument("--dna", action='store_true', dest='with_dna', help="DNA is present in the system")
    # parser.add_argument("-d", "--outdir", help="Output directory")

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug("Starting contact analysis with MDAnalysis.")

    # Check that group1 and group2 are provided
    if group1 is None or group2 is None:
        logger.error("Both group1 and group2 must be specified for interaction analysis.")
        sys.exit(1)

    # Parameters
    cutoff = 5.0  # Distance cutoff in Å

    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[0]

    # Define selections
    selections = []
    for group in [group1, group2]:
        if re.match(r"^[A-Za-z]$", group):
            logger.debug(f"MDA - contact analysis - Group {group} identified as chain.")
            selections.append(u.select_atoms(f"segid {group}"))
        elif group == "DNA" and dna:
            logger.debug(f"MDA - contact analysis - Group {group} identified as DNA.")
            selections.append(u.select_atoms("nucleic")) 
        elif group == "Protein":
            logger.debug(f"MDA - contact analysis - Group {group} identified as Protein.")
            protein_sel = u.select_atoms("protein")
            # Get segids of the protein, to avoid leaving out non-standard residues like PTR
            prot_segids = [seg.segid for seg in protein_sel.segments]
            logger.debug(f"Protein segids identified: {prot_segids}")
            selections.append(u.select_atoms("segid " + " ".join(prot_segids)))

    if len(selections) != 2:
        logger.warning("Error in defining selections for contact analysis.")

    sel1, sel2 = selections

    # Pairwise distance calculation within cutoff 5.0 Å 
    if cores > 1:
        logger.debug(f"MDA - contact analysis - Running parallel Dask for contact pairs with {cores} cores")
        contacts_pairs = an_perblock(logger, tprfile, xtcfile, cores, contact_pairs, start, stop, cutoff, sel1, sel2)
    else:
        contacts_pairs = contact_pairs(u, start, stop, cutoff, sel1, sel2)
    
    contacts_df = pd.DataFrame(contacts_pairs, 
                               columns=["frame","sele1_resid","sele1_resnm","sele1_chain","sele2_resid","sele2_resnm","sele2_chain","distance"])
    contacts_df.to_csv(outfiles["pairs_full_csv"], index=False)


if __name__ == "__main__":
    main()