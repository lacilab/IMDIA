#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scripts.utils.helpers import setup_logger


def main(logfile: str, csv: Path, ndx: Path, outdir: Path):
    """
    Plot gromos clustering results.

    Arguments:
        logfile   : Path to the log file.
        csv       : Path to .csv file containing PCA results.
        ndx       : Path to GROMACS .ndx file generating during gromos clustering (e.g. analysis/clust/gromos_clusters.ndx)
        outdir    : Output directory for results (e.g. rep1/analysis/clust).
    """

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    df_pca = pd.read_csv(csv, sep=",", header='infer')

    # Parse clusters.ndx and count cluster sizes
    f_clusters = {}
    current_cluster = None
    with open(ndx, 'r') as f:
        for line in f:
            if line.startswith('[Cluster_'):
                cluster_str = line.split('Cluster_')[1].rstrip(']\n\r')
                current_cluster = int(cluster_str)
            elif line.strip() and current_cluster is not None:
                frames = [int(frame) for frame in line.split()]
                for frame in frames:
                    f_clusters[frame] = current_cluster
    logger.debug(f"Parsed {len(f_clusters)} frames, classified into {len(set(f_clusters.values()))} clusters, from {ndx}.")

    # Map time to frames and assign clusters
    df_clusters = pd.DataFrame({'Frame': list(f_clusters.keys()), 'Cluster': list(f_clusters.values())})
    tstep = df_pca['Time(ps)'].iloc[1] - df_pca['Time(ps)'].iloc[0]
    df_clusters['Time(ps)'] = (df_clusters['Frame'] - 1) * tstep  # Frames start at 1, time starts at 0

    df_pca_ext = pd.merge(df_pca, df_clusters, on='Time(ps)', how='inner')
    df_pca_ext.rename(columns={'Cluster': 'Gromos_Cluster'}, inplace=True)
    df_reordered = df_pca_ext[['Time(ps)', 'Frame', 
                               'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 
                               'KMeans_Cluster', 'Hierarchical_Cluster',  'Gromos_Cluster']]
    csv_ext = outdir / 'pca_clusters_ext.csv'
    df_reordered.to_csv(csv_ext, index=False)
    logger.info(f"Extended PCA data with Gromos clusters saved to {csv_ext}")

    # Plot PC1 vs PC2, colored by cluster
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(df_pca['PC1'], df_pca['PC2'], c=df_clusters['Cluster'], cmap='tab20', s=20)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PC1 vs PC2 Colored by Gromos Cluster (0.1 nm neighbor cutoff)')
    plt.grid(True)
    plt.legend(*scatter.legend_elements(), title="Cluster ID", loc='best')
    plt.savefig(f'{outdir}/pca_gromos_clusters.png')
    plt.close()
    logger.info(f"Plot of PCA colored by gromos clusters saved to {outdir}/pca_gromos_clusters.png")

    # Plot PC1 vs PC2, first 10 clusters in blue, others in gray
    colors = ['blue' if c <= 10 else 'gray' for c in df_clusters['Cluster']]
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(df_pca['PC1'], df_pca['PC2'], c=colors, s=20)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PC1 vs PC2 Colored by Gromos Cluster (Clusters 1-10 in Blue)')
    plt.grid(True)
    plt.savefig(f'{outdir}/pca_gromos_clusters_top10.png')
    plt.close()
    logger.info(f"Plot of PCA colored by first 10 gromos clusters saved to {outdir}/pca_gromos_clusters_top10.png")

    # Bar plot of cluster frequencies
    cluster_ids = sorted(set(df_clusters['Cluster']))
    sizes = [df_clusters[df_clusters['Cluster'] == cid].shape[0] for cid in cluster_ids]
    total_frames = sum(sizes)  # Total number of frames across all clusters
    frequencies = [ (size / total_frames) for size in sizes]  # Convert to percentage

    plt.bar(cluster_ids, frequencies)
    plt.xlabel('Cluster ID')
    plt.ylabel('Frequency (% of frames)')
    plt.title('Cluster Frequencies')
    plt.tight_layout()
    plt.grid(True, axis='y')
    plt.savefig(f'{outdir}/gromos_clusters_freq.png')
    plt.close()
    logger.info(f"Bar plot of gromos cluster frequencies saved to {outdir}/gromos_clusters_freq.png")

    # Cluster ID vs time plot
    time, cluster_id = np.loadtxt(f'{outdir}/gromos_clust_ids.xvg', comments=['@','#', "&"], unpack=True)
    plt.scatter(time, cluster_id, c=cluster_id, cmap='tab10', s=20)
    plt.ylabel('Gromos Cluster ID')
    plt.gca().set_yticks(range(min(cluster_ids), max(cluster_ids) + 1, 1))
    plt.xlabel('Time (ps)')
    plt.title('Gromos Cluster ID vs Time')
    plt.savefig(f'{outdir}/gromos_clust_ids.png')
    plt.close()
    logger.info(f"Plot of gromos cluster ID vs time saved to {outdir}/gromos_clust_ids.png")

    # RMSD distance distribution plot
    x, y = np.loadtxt(f'{outdir}/gromos_clusters_rmsd_dist.xvg', comments=['@','#'], unpack=True)
    plt.plot(x, y)
    plt.xlabel('RMSD (nm)')
    plt.ylabel('Frame pair counts')
    plt.title('RMSD Distance Distribution')
    plt.grid(True)
    plt.savefig(f'{outdir}/gromos_rmsd_dist.png')
    plt.close()
    logger.info(f"RMSD distance distribution plot saved to {outdir}/gromos_rmsd_dist.png")



if __name__ == "__main__":
    main()