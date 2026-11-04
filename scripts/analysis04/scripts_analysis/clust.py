import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score
from scripts.utils.helpers import setup_logger


def main(logfile: Path, outdir: Path, csv: Path):
    """
    Run clustering with k-means and hierarchical clustering algorithms in scikit-learn and plot results.

    Arguments:
        logfile   : Path to the log file.
        outdir    : Output directory for results (e.g. rep1/analysis/clust).
        csv       : Path to .csv file containing PCA results.
    """

    if not logfile.is_file():
      print(f"Log file {logfile} does not exist.")
      sys.exit(1)
    
    logger = setup_logger("analysis", logfile)

    if not csv.is_file():
      logger.error(f"CSV file {csv} passed for clustering does not exist.")
      sys.exit(1)

    outdir.mkdir(parents=True, exist_ok=True)

    # Expected output files
    filenames = [
        'pca_clusters.csv',
        'clust_kmeans.png',
        'clust_kmeans_sizes.png',
        'clust_kmeans_time.png',
        'clust_hierarchical.png',
        'clust_hierarchical_sizes.png',
        'clust_hierarchical_time.png',
        'clust_silhouette.png',
        'summary.txt'
    ]
    outdir = Path(outdir)
    outfiles = [outdir / n for n in filenames]
    
    # Plot first two PCs
    df = pd.read_csv(csv)
    pcs = df[['PC1', 'PC2']].values
    total_frames = pcs.shape[0]
    logger.debug(f"Loaded {total_frames} frames, PC1-PC2 space")
    
    n_clusters_range = range(2, 8)
    
    # KMeans
    kmeans_results = {}
    for n in n_clusters_range:
        kmeans = KMeans(n_clusters=n, n_init=20).fit(pcs)
        labels = kmeans.labels_
        score = silhouette_score(pcs, labels)
        kmeans_results[n] = {'labels': labels, 'score': score, 'center': kmeans.cluster_centers_}
        logger.debug(f"KMeans k={n}: silhouette={score:.3f}")
    
    # Hierarchical (Ward linkage)
    hclust_results = {}
    for n in n_clusters_range:
        hclust = AgglomerativeClustering(n_clusters=n, linkage='ward').fit(pcs)
        labels = hclust.labels_
        score = silhouette_score(pcs, labels)
        hclust_results[n] = {'labels': labels, 'score': score}
        logger.debug(f"Hierarchical k={n}: silhouette={score:.3f}")
    
    best_kmeans = max(kmeans_results.keys(), key=lambda k: kmeans_results[k]['score'])
    best_hclust = max(hclust_results.keys(), key=lambda k: hclust_results[k]['score'])
    
    logger.debug(f"Best KMeans: k={best_kmeans}, score={kmeans_results[best_kmeans]['score']:.3f}")
    logger.debug(f"Best Hierarchical: k={best_hclust}, score={hclust_results[best_hclust]['score']:.3f}")

    # Save cluster labels to CSV
    df['KMeans_Cluster'] = kmeans_results[best_kmeans]['labels']
    df['Hierarchical_Cluster'] = hclust_results[best_hclust]['labels']
    df.to_csv(outdir / 'pca_clusters.csv', index=False)
    logger.debug(f"Saved PCA coordinates with cluster labels to '{outdir / 'pca_clusters.csv'}'")
    
    #### Plotting   
    # KMeans
    km_labels = kmeans_results[best_kmeans]['labels']
    plt.scatter(pcs[:,0], pcs[:,1], c=km_labels, cmap='tab10', alpha=0.7)
    plt.scatter(kmeans_results[best_kmeans]['center'][:,0], 
                     kmeans_results[best_kmeans]['center'][:,1], 
                     c='red', marker='*', s=200, ec='black', linewidth=1)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title(f'KMeans k={best_kmeans} (Centroids)')
    plt.savefig(f"{outdir}/clust_kmeans.png")
    plt.close()

    # Cluster sizes (KMeans)
    unique, counts = np.unique(km_labels, return_counts=True)
    plt.bar(unique, counts)
    plt.xlabel('Cluster ID')
    plt.ylabel('Number of Frames')
    plt.title('KMeans Cluster Sizes')
    plt.savefig(f"{outdir}/clust_kmeans_sizes.png")
    plt.close()

    # Cluster ID vs time for KMeans
    time_frames = np.arange(1, len(km_labels) + 1)
    plt.scatter(time_frames, km_labels, c=km_labels, cmap='tab10', alpha=0.7, s=20)
    plt.xlabel('Frame number')
    plt.ylabel('KMeans Cluster ID')
    plt.gca().set_yticks(np.arange(min(km_labels), max(km_labels)+1, 1))
    plt.title(f'KMeans k={best_kmeans} vs Time')
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{outdir}/clust_kmeans_time.png")
    plt.close()

    # Hierarchical  
    hc_labels = hclust_results[best_hclust]['labels']
    unique_labels = np.unique(hc_labels)
    medoids = []

    for label in unique_labels:
        cluster_points = pcs[hc_labels == label]
        # Find point closest to cluster mean (medoid)
        centroid = cluster_points.mean(axis=0)
        distances = np.linalg.norm(cluster_points - centroid, axis=1)
        medoid_idx = np.argmin(distances)
        medoids.append(cluster_points[medoid_idx])

    medoids = np.array(medoids)

    plt.scatter(pcs[:,0], pcs[:,1], c=hc_labels, cmap='tab10', alpha=0.7)
    plt.scatter(medoids[:,0], medoids[:,1], c='red', marker='*', s=200, ec='black', linewidth=1)
    plt.title(f'Hierarchical k={best_hclust} (Medoids)')
    plt.savefig(f"{outdir}/clust_hierarchical.png")
    plt.close()

    # Cluster sizes (Hierarchical)
    clust_id, counts = np.unique(hc_labels, return_counts=True)
    plt.bar(clust_id, counts)
    plt.xlabel('Cluster ID')
    plt.ylabel('Number of Frames')
    plt.title('Hierarchical Cluster Sizes')
    plt.savefig(f"{outdir}/clust_hierarchical_sizes.png")
    plt.close()

    # Cluster ID vs time for Hierarchical 
    plt.scatter(time_frames, hc_labels, c=hc_labels, cmap='tab10', alpha=0.7, s=20)
    plt.xlabel('Frame number')
    plt.ylabel('Hierarchical Cluster ID')
    plt.gca().set_yticks(np.arange(min(hc_labels), max(hc_labels)+1, 1))
    plt.title(f'Hierarchical k={best_hclust} (Ward) vs Time')
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{outdir}/clust_hierarchical_time.png")
    plt.close()
    
    # Silhouette scores
    x = np.arange(2, 8)
    plt.plot(x, [kmeans_results[n]['score'] for n in x], 'o-', label='KMeans', linewidth=2)
    plt.plot(x, [hclust_results[n]['score'] for n in x], 's-', label='Hierarchical', linewidth=2)
    plt.xlabel('k')
    plt.ylabel('Silhouette Score')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.title('Cluster Quality')
    plt.savefig(f"{outdir}/clust_silhouette.png")
    plt.close()

    #### Text summary
    with open(outdir / 'summary.txt', 'w') as f:
        # Header
        f.write("Clustering Summary\n")
        f.write("="*50 + "\n\n")

        # Silhouette scores
        for method, results, best_k in [("KMeans", kmeans_results, best_kmeans), 
                                        ("Hierarchical (Ward)", hclust_results, best_hclust)]:
            f.write(f"{method} Silhouette Scores:\n")
            for n in sorted(results.keys()):
                f.write(f"  k={n}: {results[n]['score']:.3f}\n")
            f.write(f"Best k for {method}: {best_k} with silhouette score {results[best_k]['score']:.3f}\n\n")

        # Cluster statistics for both methods
        for method, results, best_k in [("KMeans", kmeans_results, best_kmeans), 
                                        ("Hierarchical (Ward)", hclust_results, best_hclust)]:
            f.write(f"Cluster statistics ({method}):\n")
            f.write(f"Best k: {best_k}\n")
            
            labels = results[best_k]['labels']
            unique_labels = np.unique(labels)
            
            for label in sorted(unique_labels):
                cluster_mask = labels == label
                cluster_frames = np.sum(cluster_mask)
                cluster_pct = cluster_frames / total_frames * 100
                
                # Medoid (closest to centroid) as representative
                cluster_pcs = pcs[cluster_mask]
                centroid = cluster_pcs.mean(axis=0)
                distances = np.linalg.norm(cluster_pcs - centroid, axis=1)
                rep_idx = np.argmin(distances)
                rep_frame = np.where(cluster_mask)[0][rep_idx] + 1  # 1-based frame indexing
                
                f.write(f"Cluster {label}: {cluster_frames} frames ({cluster_pct:.1f}%)\n")
                f.write(f"Representative frame index: {rep_frame}\n")
            
            f.write("\n")

    # Check all expected outputs were generated
    missing = [f.name for f in outfiles if not f.is_file()]
    if missing:
        logger.warning("K-means and hierarchical clustering did not complete successfully. Some expected outputs are missing.")
        logger.warning(f"Missing {len(missing)} file(s): {', '.join(missing)}")
        logger.warning(f"Check log {logfile}")
    else:
        logger.info(f"K-means and hierarchical clustering completed successfully. Output files are in {outdir}.")
    

if __name__ == "__main__":
    main()
