import pickle
import sys
from MDAnalysis.analysis import rms, align
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# Arguments
datadir = sys.argv[1]
universe = sys.argv[2]
rmsd = sys.argv[3]
rmsd_cats = sys.argv[4]
rmsd_bb = sys.argv[5]
radius_axis = sys.argv[6]
radius_rmsd = sys.argv[7]
rmsf_prot_domains = sys.argv[8]
rmsf_dna_res = sys.argv[9]
rmsd_csv = sys.argv[10]
radius_axis_csv = sys.argv[11]
rmsf_prot_csv = sys.argv[12]
rmsf_dna_csv = sys.argv[13]

os.makedirs(datadir, exist_ok=True)

with open(universe, 'rb') as f:
    u = pickle.load(f)
u.trajectory[0]

# Select atoms by chain
dimer = u.select_atoms("segid A B")  
monoA = u.select_atoms("segid A")   
monoB = u.select_atoms("segid B")   
dna = u.select_atoms("segid C D")   
dnaC = u.select_atoms("segid C") 
dnaD = u.select_atoms("segid D") 

########## RMSD ##########
rmsd_analysis = rms.RMSD(u, select='backbone', groupselections=['name CA', 'protein'])
rmsd_analysis.run()
print(rmsd_analysis.results.rmsd.shape)
rmsd_df = pd.DataFrame(rmsd_analysis.results.rmsd[:, 2:],
                columns=['Backbone', 'C-alphas', 'Protein'],
                index=rmsd_analysis.results.rmsd[:, 1])
rmsd_df.index.name = 'Time (ps)'

rmsd_df.to_csv(rmsd_csv)

rmsd_df.plot(title='RMSD')
plt.ylabel('RMSD (Å)')
plt.savefig(rmsd)
plt.close()

rmsd_df_bb = pd.DataFrame(rmsd_analysis.results.rmsd[:, 2],
                columns=['Backbone'],
                index=rmsd_analysis.results.rmsd[:, 1])
rmsd_df_bb.index.name = 'Time (ps)'

rmsd_df_ca = pd.DataFrame(rmsd_analysis.results.rmsd[:, 3],
                    columns=['C-alphas'],
                    index=rmsd_analysis.results.rmsd[:, 1])
rmsd_df_ca.index.name = 'Time (ps)'

rmsd_df_prot = pd.DataFrame(rmsd_analysis.results.rmsd[:, 4],
                    columns=['Protein'],
                    index=rmsd_analysis.results.rmsd[:, 1])
rmsd_df_prot.index.name = 'Time (ps)'

fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # 1 row, 3 columns
rmsd_df_bb.plot(title='RMSD Backbone', ax=axes[0], legend=False)
rmsd_df_ca.plot(title='RMSD C-alphas', ax=axes[1], legend=False, c='orange')
rmsd_df_prot.plot(title='RMSD Protein', ax=axes[2], legend=False, c='green')
plt.tight_layout()
plt.savefig(rmsd_cats)
plt.close()

rmsd_df_bb.plot(title='RMSD Backbone', legend=False)
plt.ylabel('RMSD (Å)')
plt.savefig(rmsd_bb)
plt.close()

########## Radius of gyration ##########
def radgyr(atomgroup, masses, total_mass=None):
    # coordinates change for each frame
    coordinates = atomgroup.positions
    center_of_mass = atomgroup.center_of_mass()

    # get squared distance from center
    ri_sq = (coordinates-center_of_mass)**2
    # sum the unweighted positions
    sq = np.sum(ri_sq, axis=1)
    sq_x = np.sum(ri_sq[:,[1,2]], axis=1) # sum over y and z
    sq_y = np.sum(ri_sq[:,[0,2]], axis=1) # sum over x and z
    sq_z = np.sum(ri_sq[:,[0,1]], axis=1) # sum over x and y

    # make into array
    sq_rs = np.array([sq, sq_x, sq_y, sq_z])

    # weight positions
    rog_sq = np.sum(masses*sq_rs, axis=1)/total_mass
    # square root and return
    return np.sqrt(rog_sq)

result = []
for frame in u.trajectory:
    result.append(radgyr(atomgroup=dimer,
                        masses=dimer.masses,
                        total_mass=np.sum(dimer.masses)))

result = np.asarray(result).T
labels = ['all', 'x-axis', 'y-axis', 'z-axis']

# Save as CSV
df = pd.DataFrame(result, index=labels, columns=np.arange(len(u.trajectory))).T
df.index.name = "Frame"
df.to_csv(radius_axis_csv)

for col, label in zip(result, labels):
    plt.plot(col, label=label)
plt.legend()
plt.ylabel('Radius of gyration (Å)')
plt.xlabel('Frame')
plt.savefig(radius_axis)
plt.close()

plt.figure(figsize=(8, 5))
plt.scatter(rmsd_analysis.results.rmsd[:, 2], result[0,:],s=10, alpha=0.4, label="RMSD vs Radius of gyration")
plt.xlabel("RMSD Backbone (Å)")
plt.ylabel("Radius of gyration (Å)")
plt.title("RMSD vs Radius of gyration")
plt.legend()
plt.savefig(radius_rmsd)
plt.close()

########## RMSF ##########
average = align.AverageStructure(u, u, select='protein and name CA',
                            ref_frame=0).run()
ref = average.results.universe

aligner = align.AlignTraj(u, ref,
                    select='protein and name CA',
                    in_memory=True).run()

monoA_ca = monoA.select_atoms('name CA')
RA = rms.RMSF(monoA_ca).run()
monoB_ca = monoB.select_atoms('name CA')
RB = rms.RMSF(monoB_ca).run()

# Save Protein RMSF to CSV
rmsf_prot = pd.DataFrame({
    'Residue': monoA_ca.resids,
    'RMSF_A': RA.results.rmsf,
    'RMSF_B': RB.results.rmsf
})
rmsf_prot.to_csv(rmsf_prot_csv, index=False)

# Plot
plt.plot(monoA_ca.resids + 135, RA.results.rmsf, label='RMSF monomer A')
plt.plot(monoB_ca.resids + 135, RB.results.rmsf, label='RMSF monomer B')
plt.axvspan(137, 317, zorder=0, alpha=0.2, color='blue', label='CCD')
plt.axvspan(318, 488, zorder=0, alpha=0.2, color='red', label='DBD')
plt.axvspan(489, 576, zorder=0, alpha=0.2, color='green', label='LD')
plt.axvspan(577, 683, zorder=0, alpha=0.2, color='yellow', label='SH2')
plt.axvspan(684, 708, zorder=0, alpha=0.2, color='pink', label='Tail')
plt.xlabel('Residue number (adjusted)')
plt.ylabel('RMSF (Å)')
plt.legend()
plt.savefig(rmsf_prot_domains)
plt.close()

dnaC_p = dnaC.select_atoms('name P')
RC = rms.RMSF(dnaC_p).run()
dnaD_p = dnaD.select_atoms('name P')
RD = rms.RMSF(dnaD_p).run()

# Save DNA RMSF to CSV
rmsf_dna = pd.DataFrame({
    'Residue': dnaC_p.resids,
    'RMSF_C': RC.results.rmsf,
    'RMSF_D': RD.results.rmsf
})
rmsf_dna.to_csv(rmsf_dna_csv, index=False)

# Plot
plt.plot(dnaC_p.resids, RC.results.rmsf, label='RMSF DNA C')
plt.plot(dnaD_p.resids, RD.results.rmsf, label='RMSF DNA D')
plt.xlabel('Residue number')
plt.ylabel('RMSF (Å)')
plt.legend()
plt.savefig(rmsf_dna_res)
plt.close()