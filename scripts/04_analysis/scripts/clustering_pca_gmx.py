#!/usr/bin/env python3

import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd

# Arguments
parser = argparse.ArgumentParser(description="Project sklearn clustering results on GROMACS PCA data.")
parser.add_argument("-s", "--timescale", type=str, help="Path to the 2D projection timescale of the PCA (CSV file).")
parser.add_argument("-l","--labels", type=str, help="Clustering labels output by sklearn. Expected format: algorithm_labels.csv")
parser.add_argument("-d", "--outdir", type=str, help="Output directory")

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(0)

args = parser.parse_args()
pca_data_file = args.timescale
clustering_labels_file = args.labels
outputdir = args.outdir

# TODO check files and directories exist
# TODO check file extensions are correct

# Extract clustering algorithm from filename
clustering_alg = clustering_labels_file.split('_')[0]

# Load csv files
pca_data = pd.read_csv(pca_data_file, sep=",", header=None, skipfooter=1, names=['PCA1', 'PCA2', 'Time'])
clustering_labels = pd.read_csv(clustering_labels_file)
merged_data = pd.concat([pca_data, clustering_labels['Cluster']], axis=1)
merged_data.to_csv(f'{outputdir}/pca_clusters_{clustering_alg}.csv', index=False) 

# data = pd.read_csv(f'{outputdir}/pca_kmeans.csv')
plt.figure(figsize=(10, 6))
scatter = plt.scatter(
    merged_data['PCA1'], merged_data['PCA2'], # convert from nm to Å
    c=merged_data['Cluster'], cmap='tab20', alpha=0.6
)

plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=10)) 
# plt.colorbar(scatter, label='Cluster')
plt.title(f'GROMACS PCA Projection with {clustering_alg} Clusters')
plt.xlabel('PC1 (Å)')
plt.ylabel('PC2 (Å)')
plt.savefig(f'{outputdir}/pca_gmx_clusters_{clustering_alg}.jpg', dpi=150)
plt.close()