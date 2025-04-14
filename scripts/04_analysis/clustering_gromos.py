import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import sys

# Arguments
pca_data_file = sys.argv[1]
clusters_index_file = sys.argv[2]
outputdir = sys.argv[3]

pca_data = pd.read_csv(pca_data_file, sep=r"\s+", header=None, skipfooter=1, names=['PC1', 'PC2', 'Time'])

# Parse clusters.ndx and count cluster sizes
frame_to_cluster = {}
cluster_sizes = {}
current_cluster = None
with open(clusters_index_file, 'r') as f:
    for line in f:
        if line.startswith('[Cluster_'):
            current_cluster = int(line.strip()[9:-1])  # Extract cluster number
            cluster_sizes[current_cluster] = 0 
        elif line.strip() and current_cluster is not None:
            frames = [int(frame) for frame in line.split()]
            for frame in frames:
                frame_to_cluster[frame] = current_cluster
                cluster_sizes[current_cluster] += 1 

# Map time to frames and assign clusters
clusters = []
for t in pca_data['Time']:
    frame = int(t/100)
    cluster_id = frame_to_cluster.get(frame, 0) 
    clusters.append(cluster_id)

# Plot PC1 vs PC2, colored by cluster
plt.figure(figsize=(8, 6))
scatter = plt.scatter(pca_data['PC1'], pca_data['PC2'], c=clusters, cmap='tab20', s=20)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PC1 vs PC2 Colored by Gromos Cluster')
plt.grid(True)
plt.savefig(f'{outputdir}/pca_gromos_clusters.png')
plt.close()

# Plot PC1 vs PC2, first 10 clusters in blue, others in gra
colors = ['blue' if c <= 10 else 'gray' for c in clusters]
plt.figure(figsize=(8, 6))
scatter = plt.scatter(pca_data['PC1'], pca_data['PC2'], c=colors, s=20)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PC1 vs PC2 Colored by Gromos Cluster (1-10 in Blue)')
plt.grid(True)
plt.savefig(f'{outputdir}/pca_gromos_clusters_top10.png')
plt.close()

# Bar plot of cluster frequencies
cluster_ids = sorted(cluster_sizes.keys())
sizes = [cluster_sizes[cid] for cid in cluster_ids]
total_frames = sum(sizes)  # Total number of frames across all clusters
frequencies = [ (size / total_frames) for size in sizes]  # Convert to percentage

plt.bar(cluster_ids, frequencies, color='skyblue')
plt.xlabel('Cluster ID')
plt.ylabel('Frequency')
plt.title('Cluster Frequencies')
plt.grid(True, axis='y')
plt.tight_layout()
plt.savefig(f'{outputdir}/gromos_clusters_freq.png')
plt.close()

# Bar plot of cumulative cluster frequencies
cumulative_frequencies = [sum(frequencies[:i+1]) for i in range(len(frequencies))]  # Cumulative sum

plt.bar(cluster_ids, cumulative_frequencies, color='skyblue')
plt.xlabel('Cluster ID')
plt.ylabel('Cumulative Frequency (% of Total Frames)')
plt.title('Cumulative Cluster Frequencies')
plt.grid(True, axis='y')
plt.tight_layout()
plt.savefig(f'{outputdir}/gromos_clusters_freq_cum.png')
plt.close()