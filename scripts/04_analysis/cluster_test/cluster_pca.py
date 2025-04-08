import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd

CONFIG = {
    'output_dir': 'umap_analysis_results',
}

pca_data = pd.read_csv('2dproj_timescale.csv', sep=r"\s+", header=None, skipfooter=1, names=['PCA1', 'PCA2', 'Time'])
clustering_labels = pd.read_csv(f"{CONFIG['output_dir']}/kmeans_labels.csv")
merged_data = pd.concat([pca_data, clustering_labels['Cluster']], axis=1)
merged_data.to_csv(f'{CONFIG['output_dir']}/pca_with_clusters.csv', index=False) 

data = pd.read_csv(f'{CONFIG['output_dir']}/pca_kmeans.csv')
plt.figure(figsize=(10, 6))
scatter = plt.scatter(
    data['PCA1'], data['PCA2'],
    c=data['Cluster'], cmap='tab20', alpha=0.6
)

plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=10))  # Adjust `nbins` to control the number of ticks
plt.colorbar(scatter, label='Cluster')
plt.title('PCA Projection with Clusters')
plt.xlabel('PC 1')
plt.ylabel('PC 2')
plt.savefig('pca_kmeans_overlay.png', dpi=300)
plt.close()