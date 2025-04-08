import sys
import os
import pickle
import MDAnalysis as mda
from MDAnalysis.analysis import pca, align
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

# Arguments
datadir = sys.argv[1]
universe = sys.argv[2]
pc_cumulative = sys.argv[3]
pc_pairgrid = sys.argv[4]
pc_cosine = sys.argv[5]

os.makedirs(datadir, exist_ok=True)

with open(universe, 'rb') as f:
    u = pickle.load(f)
u.trajectory[0]

dimer = u.select_atoms("segid A B") 
monoA = u.select_atoms("segid A")    
monoB = u.select_atoms("segid B")    
dna = u.select_atoms("segid C D")  
dnaC = u.select_atoms("segid C") 
dnaD = u.select_atoms("segid D")   


aligner = align.AlignTraj(u, u, select='protein and name CA',
                          in_memory=True).run()

pc = pca.PCA(u, select='protein and name CA',
             align=True, mean=None,
             n_components=None).run()

backbone = u.select_atoms('protein and name CA')
n_bb = len(backbone)
print('There are {} backbone atoms in the analysis'.format(n_bb))
print(pc.results.p_components.shape)

# Plot cumulative variance captured by PCs
plt.plot(pc.results.cumulated_variance[:15])
plt.xlabel('Principal component')
plt.ylabel('Cumulative variance')
plt.savefig(pc_cumulative)
plt.close()

# Plot pairgrid
transformed = pc.transform(backbone, n_components=3)
df_pca = pd.DataFrame(transformed,
                  columns=['PC{}'.format(i+1) for i in range(3)])
df_pca['Time (ps)'] = df_pca.index * u.trajectory.dt
g = sns.PairGrid(df_pca, hue='Time (ps)',
                 palette=sns.color_palette('Oranges_d',
                                           n_colors=len(df_pca)))
g.map(plt.scatter, marker='.')
plt.savefig(pc_pairgrid)
plt.close()

# Cosine content
# 0 to 1, <0.7 is OK
for i in range(3):
    cc = pca.cosine_content(transformed, i)
    print(f"Cosine content for PC {i+1} = {cc:.3f}")

melted = pd.melt(df_pca, id_vars=["Time (ps)"],
                 var_name="PC",
                 value_name="Value")
g = sns.FacetGrid(melted, col="PC")
g.map(sns.lineplot,
      "Time (ps)", # x-axis
      "Value", # y-axis
      errorbar=None) # no confidence interval
plt.savefig(pc_cosine)
plt.close()
