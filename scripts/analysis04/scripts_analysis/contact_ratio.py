#!/usr/bin/env python3

import sys
import re
import pickle
from MDAnalysis.analysis.contacts import Contacts
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scripts.utils.helpers import setup_logger
from scripts.utils.parallel import an_perblock


def main(logfile : str = None, universe: str = None, 
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

    logger.debug(f"MDA - contact analysis - Native contact ratio calculation between groups: {group1} and {group2}")
    cont = Contacts(u, select=(sel1, sel2), refgroup=(sel1, sel2), radius=cutoff)
    cont.run(start=start, stop=stop, verbose=True)
    dist_results = cont.results.timeseries

    df_timeseries = pd.DataFrame(dist_results)
    df_timeseries.to_csv(outfiles["contacts_timeseries_csv"], index=False)

    average_contacts = np.mean(dist_results[:, 1])

    # Timeseries plot
    plt.figure(figsize=(10, 6))
    plt.plot(dist_results[:, 0], dist_results[:, 1], label='Contacts', color='blue')
    plt.axhline(y=average_contacts, color='red', linestyle='--', label='Average Contacts')
    plt.xlabel('Frame')
    plt.ylabel('Fraction of Contacts')
    plt.title('Native Contacts, average = {:.2f}'.format(average_contacts))
    plt.legend()
    plt.savefig(outfiles["contacts_timeseries_png"])
    plt.close()

if __name__ == "__main__":
    main()