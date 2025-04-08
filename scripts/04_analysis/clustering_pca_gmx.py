import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import sys

# Arguments
pca_data_file = sys.argv[1]
clustering_labels_file = sys.argv[2]
clustering_alg = sys.argv[3]
outputdir = sys.argv[4]

pca_data = pd.read_csv(pca_data_file, sep=r"\s+", header=None, skipfooter=1, names=['PCA1', 'PCA2', 'Time'])
clustering_labels = pd.read_csv(clustering_labels_file)
merged_data = pd.concat([pca_data, clustering_labels['Cluster']], axis=1)
merged_data.to_csv(f'{outputdir}/pca_clusters_{clustering_alg}.csv', index=False) 

# data = pd.read_csv(f'{outputdir}/pca_kmeans.csv')
plt.figure(figsize=(10, 6))
scatter = plt.scatter(
    merged_data['PCA1'], merged_data['PCA2'],
    c=merged_data['Cluster'], cmap='tab20', alpha=0.6
)

plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=10)) 
plt.colorbar(scatter, label='Cluster')
plt.title('PCA Projection with Clusters')
plt.xlabel('PC 1')
plt.ylabel('PC 2')
plt.savefig('pca_clusters_{clustering_alg}.png', dpi=300)
plt.close()