#!/usr/bin/env python3

import sys
import re
import pickle
import argparse
import MDAnalysis as mda
from MDAnalysis.analysis.contacts import Contacts
from MDAnalysis.lib.distances import distance_array
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scripts.utils.helpers import setup_logger


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
            contact = {
                "frame": frame,
                "sele1_resid": f"{atom1.resid}", 
                "sele1_resnm": f"{atom1.resname}", 
                "sele1_chain": f"{atom1.segid}",
                "sele2_resid": f"{atom2.resid}", 
                "sele2_resnm": f"{atom2.resname}",
                "sele2_chain": f"{atom2.segid}",
                "distance": dist_matrix[i, j]
            }
            contact_data.append(contact)
    return pd.DataFrame(contact_data)


##### MAIN FUNCTION #####

def main(logfile : str = None, universe: str = None, group1 : str = None, group2 : str = None,
         dna : bool = False, outfiles : list = []):
    # Arguments
    # parser = argparse.ArgumentParser(description="Perform contact analysis with MDAnalysis.")
    # parser.add_argument("-u", "--universe", help="Path to the Universe object")
    # parser.add_argument("--dna", action='store_true', dest='with_dna', help="DNA is present in the system")
    # parser.add_argument("-d", "--outdir", help="Output directory")
# 
    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit(0)
# 
    # args = parser.parse_args()
    # universe = args.universe
    # outdir = args.outdir

    # Output file paths
    # timeseries_full_csv = f"{outdir}/timeseries_full.csv"
    # timeseries_full_png = f"{outdir}/timeseries_full.png"
    # timeseries_dbd_csv = f"{outdir}/timeseries_dbd.csv"
    # timeseries_dbd_png = f"{outdir}/timeseries_dbd.png"
    # pairs_dbd_csv = f"{outdir}/pairs_dbd.csv"
    # pairs_full_csv = f"{outdir}/pairs_full.csv"

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug("Starting contact analysis with MDAnalysis.")

    # Check that group1 and group2 are provided
    if group1 is None or group2 is None:
        logger.error("Both group1 and group2 must be specified for interaction analysis.")
        sys.exit(1)

    # Parameters
    cutoff = 5.0  # Distance cutoff in Å
    start_frame = 1000
    stop_frame = 2001

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

    ##### Distance-based contact analysis
    # Full protein
    cont = Contacts(u, select=(sel1, sel2), refgroup=(sel1, sel2), radius=cutoff)
    cont.run(start=start_frame, stop=stop_frame, verbose=True)

    df_timeseries = pd.DataFrame(cont.results.timeseries)
    df_timeseries.to_csv(outfiles["timeseries_full_csv"], index=False)

    average_contacts = np.mean(cont.results.timeseries[:, 1])

    # Timeseries plot
    plt.figure(figsize=(10, 6))
    plt.plot(cont.results.timeseries[:, 0], cont.results.timeseries[:, 1], label='Contacts', color='blue')
    plt.axhline(y=average_contacts, color='red', linestyle='--', label='Average Contacts')
    plt.xlabel('Frame')
    plt.ylabel('Fraction of Contacts')
    plt.title('Native Contacts, average = {:.2f}'.format(average_contacts))
    plt.legend()
    plt.savefig(outfiles["timeseries_full_png"])
    plt.close()

    # DBD
    # if with_dna:
    #     cont = Contacts(u, select=(sel1_dbd, sel2), refgroup=(sel1_dbd, sel2), radius=cutoff)
    #     cont.run(start=start_frame, stop=stop_frame, verbose=True)
 
    #     df_timeseries = pd.DataFrame(cont.results.timeseries)
    #     df_timeseries.to_csv(timeseries_dbd_csv, index=False)
 
    #     average_contacts = np.mean(cont.results.timeseries[:, 1])
 
    #     # Timeseries plot
    #     plt.figure(figsize=(10, 6))
    #     plt.plot(cont.results.timeseries[:, 0], cont.results.timeseries[:, 1], label='Contacts', color='blue')
    #     plt.axhline(y=average_contacts, color='red', linestyle='--', label='Average Contacts')
    #     plt.xlabel('Frame')
    #     plt.ylabel('Fraction of Contacts')
    #     plt.title('Native Contacts, average = {:.2f}'.format(average_contacts))
    #     plt.legend()
    #     plt.savefig(timeseries_dbd_png)
    #     plt.close()


    # Full protein
    # Create DataFrame
    contacts_df = contact_pairs(u, start_frame, stop_frame, cutoff, sel1, sel2)
    contacts_df.to_csv(outfiles["pairs_full_csv"], index=False)

    # DBD
    # if with_dna:
    #     contacts_df = contact_pairs(sel1_dbd, sel2)
    #     contacts_df.to_csv(outfiles["pairs_dbd_csv"], index=False)

if __name__ == "__main__":
    main()