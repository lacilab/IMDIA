import sys
import os
import pickle
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import numpy as np
import pandas as pd

# Arguments
datadir = sys.argv[1]
universe = sys.argv[2]
tprfile = sys.argv[3]
xtcfile = sys.argv[4]
hbonds_mda_parsed = sys.argv[5]

os.makedirs(datadir, exist_ok=True)

with open(universe, 'rb') as f:
    u = pickle.load(f)
u.trajectory[0]

u_tpr = mda.Universe(tprfile,xtcfile)

# Define groups
dimer = u.select_atoms("segid A B") 
dna = u.select_atoms("segid C D")

print("Dimer atoms:", len(dimer))
print("DNA atoms:", len(dna))

# Convert AtomGroup indices to selection strings
dimer_indices = " ".join(map(str, dimer.indices))
dna_indices = " ".join(map(str, dna.indices))

# Create selection strings
dimer_sel = f"index {dimer_indices}"
dna_sel = f"index {dna_indices}"

# Define hbond analysis parameters
hbonds_dimer_to_dna = HBA(
    universe=u_tpr,
    between=[
        [dimer_sel, dna_sel], # for protein-water hbonds
    ],
    d_h_a_angle_cutoff=120,
    d_a_cutoff=3.5
)
# Run analysis, starting at frame 1000 until 2000 (end)
hbonds_dimer_to_dna.run(start=1000, stop=2001, verbose=True)
# Extract results
hbonds_array = hbonds_dimer_to_dna.results.hbonds

# Name columns and save raw data
df = pd.DataFrame(hbonds_array, columns=["frame", "donor", "hydrogen", "acceptor", "distance", "angle"])
# Adjust to 1-index
df['donor'] = (df['donor'] + 1).astype('int')
df['hydrogen'] = (df['hydrogen'] + 1).astype('int')
df['acceptor'] = (df['acceptor'] + 1).astype('int')

# Save
df.to_csv(hbonds_mda_parsed, index=False)