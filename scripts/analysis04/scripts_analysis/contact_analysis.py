#!/usr/bin/env python3

import sys
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

# Arguments
parser = argparse.ArgumentParser(description="Perform contact analysis with MDAnalysis.")
parser.add_argument("-u", "--universe", help="Path to the Universe object")
parser.add_argument("--dna", action='store_true', dest='with_dna', help="DNA is present in the system")
parser.add_argument("-d", "--outdir", help="Output directory")

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(0)

args = parser.parse_args()
universe = args.universe
outdir = args.outdir

# Output file paths
timeseries_full_csv = f"{outdir}/timeseries_full.csv"
timeseries_full_png = f"{outdir}/timeseries_full.png"
timeseries_dbd_csv = f"{outdir}/timeseries_dbd.csv"
timeseries_dbd_png = f"{outdir}/timeseries_dbd.png"
pairs_dbd_csv = f"{outdir}/pairs_dbd.csv"
pairs_full_csv = f"{outdir}/pairs_full.csv"

# Parameters
cutoff = 5.0  # Distance cutoff in Å
start_frame = 1000
stop_frame = 2001

with open(universe, 'rb') as f:
    u = pickle.load(f)
u.trajectory[0]

# Group selection
if args.with_dna:
    sel1_full = u.select_atoms("protein or resname PTR")
    sel1_dbd = u.select_atoms("protein and resid 182:353")
    sel2 = u.select_atoms("nucleic")
else:
    sel1_full = u.select_atoms("segid A")
    #sel1_dbd = u.select_atoms("protein and resid 182:353")
    sel2 = u.select_atoms("segid B")

##### Distance-based contact analysis
# Full protein
cont = Contacts(u, select=(sel1_full, sel2), refgroup=(sel1_full, sel2), radius=cutoff)
cont.run(start=start_frame, stop=stop_frame, verbose=True)

df_timeseries = pd.DataFrame(cont.results.timeseries)
df_timeseries.to_csv(timeseries_full_csv, index=False)

average_contacts = np.mean(cont.results.timeseries[:, 1])

# Timeseries plot
plt.figure(figsize=(10, 6))
plt.plot(cont.results.timeseries[:, 0], cont.results.timeseries[:, 1], label='Contacts', color='blue')
plt.axhline(y=average_contacts, color='red', linestyle='--', label='Average Contacts')
plt.xlabel('Frame')
plt.ylabel('Fraction of Contacts')
plt.title('Native Contacts, average = {:.2f}'.format(average_contacts))
plt.legend()
plt.savefig(timeseries_full_png)
plt.close()

# DBD
if args.with_dna:
    cont = Contacts(u, select=(sel1_dbd, sel2), refgroup=(sel1_dbd, sel2), radius=cutoff)
    cont.run(start=start_frame, stop=stop_frame, verbose=True)

    df_timeseries = pd.DataFrame(cont.results.timeseries)
    df_timeseries.to_csv(timeseries_dbd_csv, index=False)

    average_contacts = np.mean(cont.results.timeseries[:, 1])

    # Timeseries plot
    plt.figure(figsize=(10, 6))
    plt.plot(cont.results.timeseries[:, 0], cont.results.timeseries[:, 1], label='Contacts', color='blue')
    plt.axhline(y=average_contacts, color='red', linestyle='--', label='Average Contacts')
    plt.xlabel('Frame')
    plt.ylabel('Fraction of Contacts')
    plt.title('Native Contacts, average = {:.2f}'.format(average_contacts))
    plt.legend()
    plt.savefig(timeseries_dbd_png)
    plt.close()


##### Distance by residue pairs

def contact_pairs(sel1, sel2):
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


# Full protein
# Create DataFrame
contacts_df = contact_pairs(sel1_full, sel2)
contacts_df.to_csv(pairs_full_csv, index=False)


# DBD
if args.with_dna:
    contacts_df = contact_pairs(sel1_dbd, sel2)
    contacts_df.to_csv(pairs_dbd_csv, index=False)