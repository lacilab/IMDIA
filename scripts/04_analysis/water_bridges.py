import sys
import pickle
import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import WaterBridgeAnalysis

# Arguments
universe = sys.argv[1]
wb_csv = sys.argv[2]

with open(universe, 'rb') as f:
    u = pickle.load(f)
u.trajectory[0]


# Custom forcefield atom names
class WaterBridgeAnalysis_OtherFF(WaterBridgeAnalysis):
      DEFAULT_DONORS = {"OtherFF": tuple(set(["O", "H2", "N", "NE", "NH1", "NH2", "ND2", "SG", "NE2", "ND1", 
                                              "NZ", "OG", "OG1", "OH"]))}
      DEFAULT_ACCEPTORS = {"OtherFF": tuple(set(["O", "H2", "OC1", "OC2", "OD1", "OD2", 
                                                 "SG", "OE1", "OE2", "ND1", "NE2", "SD", "OG", "OG1", "OH"]))}


wba = WaterBridgeAnalysis_OtherFF(
    universe=u,
    selection1='protein or resname PTR',  
    selection2='nucleic',  
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
df.to_csv(wb_csv, index=False)
