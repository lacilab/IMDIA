import sys
import mdtraj as md
import numpy as np
from umap.umap_ import UMAP
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # For 3D plotting
import seaborn as sns
from typing import Tuple, List, Dict
import os
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm

# Arguments
pdb = sys.argv[1]
xtc = sys.argv[2]
outputdir = sys.argv[3]

# ============= CONFIGURATION SECTION =============
CONFIG = {
    'input_files': {
        'pdb_path': pdb,
        'xtc_path': xtc,
    },

    'binding_site': {
        'residue_numbers': list(range(1, 576)),
    },

    'umap_params': {
        'n_components': 5,
        'n_neighbors_range': [5, 10, 15, 20],  # Adjusted: Balance local and global structure
        'min_dist_range': [0.01, 0.05, 0.1],  # Kept for tight packing
    },

    'clustering_params': {
        'n_clusters_range': list(range(2, 7)),
        'forced_n_clusters': 5,  # Added: Force k=5 to match initial and Ward's separation
        'hclust_linkage': 'ward'
    },

    'output_dir': outputdir,

    'stride': 1,
    'random_seed': 42
}

def load_trajectory(config: Dict) -> Tuple[md.Trajectory, np.ndarray]:
    print("Loading trajectory...")
    traj = md.load(
        config['input_files']['xtc_path'],
        top=config['input_files']['pdb_path'],
        stride=config['stride']
    )
    
    binding_site_atoms = traj.topology.select(
        f"chainid 0 1 and residue {' '.join(map(str, config['binding_site']['residue_numbers']))}"
    )
    
    if len(binding_site_atoms) == 0:
        raise ValueError("No atoms selected for binding site. Please check residue numbers.")
    
    print(f"Loaded trajectory with {traj.n_frames} frames")
    print(f"Selected {len(binding_site_atoms)} atoms for binding site analysis")
    
    return traj, binding_site_atoms

def calculate_rmsd_matrix(traj: md.Trajectory, binding_site_atoms: np.ndarray) -> np.ndarray:
    print("Calculating RMSD matrix...")
    binding_site_traj = traj.atom_slice(binding_site_atoms)
    n_frames = binding_site_traj.n_frames
    rmsd_matrix = np.zeros((n_frames, n_frames), dtype=np.float32)
    
    for i in tqdm(range(n_frames)):
        target = binding_site_traj[i:]
        reference = binding_site_traj[i]
        rmsd_values = md.rmsd(target, reference)
        rmsd_matrix[i, i:] = rmsd_values
        rmsd_matrix[i:, i] = rmsd_values
    
    return rmsd_matrix

def optimize_umap_parameters(coords: np.ndarray, config: Dict) -> Dict:
    print("Optimizing UMAP parameters...")
    os.makedirs(config['output_dir'], exist_ok=True)
    
    best_score = float('-inf')
    best_params = {}
    results = []
    
    total_combinations = len(config['umap_params']['n_neighbors_range']) * len(config['umap_params']['min_dist_range'])
    progress = tqdm(total=total_combinations)
    
    for n_neighbors in config['umap_params']['n_neighbors_range']:
        for min_dist in config['umap_params']['min_dist_range']:
            umap = UMAP(
                n_components=config['umap_params']['n_components'],
                n_neighbors=n_neighbors,
                min_dist=min_dist,
                random_state=config['random_seed']
            )
            embedded = umap.fit_transform(coords)
            
            original_dist = pdist(coords)
            embedded_dist = pdist(embedded)
            correlation = np.corrcoef(original_dist, embedded_dist)[0,1]
            
            results.append({
                'n_neighbors': n_neighbors,
                'min_dist': min_dist,
                'correlation': correlation
            })
            
            if correlation > best_score:
                best_score = correlation
                best_params = {'n_neighbors': n_neighbors, 'min_dist': min_dist}
            
            progress.update(1)
    
    progress.close()
    
    plt.figure(figsize=(10, 6))
    results_df = pd.DataFrame(results)
    pivot_table = results_df.pivot(index='n_neighbors', columns='min_dist', values='correlation')
    sns.heatmap(pivot_table, annot=True, cmap='viridis', fmt='.3f')
    plt.title('UMAP Parameter Optimization (Distance Preservation)')
    plt.savefig(f"{config['output_dir']}/umap_optimization.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    return best_params

def compute_bic(coords: np.ndarray, config: Dict) -> int:
    print("Computing BIC for cluster selection...")
    os.makedirs(config['output_dir'], exist_ok=True)
    
    bic_scores = []
    n_clusters_range = config['clustering_params']['n_clusters_range']

    for n_clusters in n_clusters_range:
        gmm = GaussianMixture(n_components=n_clusters, random_state=config['random_seed'])
        gmm.fit(coords)
        bic = gmm.bic(coords)
        bic_scores.append(bic)

    bic_scores = np.array(bic_scores)
    deltas = np.diff(bic_scores)
    second_deltas = np.diff(deltas)

    min_bic_idx = np.argmin(bic_scores)
    significant_drop = np.where(deltas < -np.mean(np.abs(deltas)))[0]
    if len(significant_drop) > 0:
        post_drop_indices = significant_drop[significant_drop < len(second_deltas)]
        if len(post_drop_indices) > 0:
            second_deriv_post_drop = second_deltas[post_drop_indices]
            optimal_idx = post_drop_indices[np.argmax(second_deriv_post_drop)] + 1
        else:
            optimal_idx = min_bic_idx
    else:
        optimal_idx = min_bic_idx

    optimal_n_clusters = n_clusters_range[optimal_idx]
    print(f"Selected optimal number of clusters: {optimal_n_clusters} (BIC: {bic_scores[optimal_idx]:.2f})")

    fig, ax1 = plt.subplots(figsize=(12, 8))
    ax1.plot(n_clusters_range, bic_scores, marker='o', label='BIC', color='blue')
    ax1.axvline(x=optimal_n_clusters, color='red', linestyle='--', label=f'Optimal k={optimal_n_clusters}')
    ax1.set_xlabel('Number of Clusters')
    ax1.set_ylabel('BIC', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.legend(loc='upper left')

    ax2 = ax1.twinx()
    ax2.plot(n_clusters_range[1:], deltas, marker='s', label='1st Derivative', color='green')
    ax2.set_ylabel('Derivatives', color='green')
    ax2.tick_params(axis='y', labelcolor='green')

    ax3 = ax1.twinx()
    ax3.spines['right'].set_position(('outward', 60))
    ax3.plot(n_clusters_range[2:], second_deltas, marker='^', label='2nd Derivative', color='orange')
    ax3.set_ylabel('2nd Derivative', color='orange')
    ax3.tick_params(axis='y', labelcolor='orange')

    fig.legend(loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)
    plt.title('BIC and Derivatives for Cluster Selection')
    plt.savefig(f"{config['output_dir']}/bic_selection_with_derivatives.png", dpi=300, bbox_inches='tight')
    plt.close()

    return optimal_n_clusters

def perform_clustering(embedded_coords: np.ndarray, config: Dict) -> Dict:
    print("Performing clustering...")
    # Check if forced_n_clusters is set; otherwise, use BIC
    if 'forced_n_clusters' in config['clustering_params'] and config['clustering_params']['forced_n_clusters'] is not None:
        optimal_n_clusters = config['clustering_params']['forced_n_clusters']
        print(f"Forcing number of clusters to: {optimal_n_clusters}")
    else:
        optimal_n_clusters = compute_bic(embedded_coords, config)

    kmeans = KMeans(
        n_clusters=optimal_n_clusters,
        random_state=config['random_seed'],
        n_init=20
    )
    kmeans_labels = kmeans.fit_predict(embedded_coords)
    kmeans_silhouette = silhouette_score(embedded_coords, kmeans_labels)
    print(f"K-means Silhouette Score: {kmeans_silhouette:.3f}")

    hclust = AgglomerativeClustering(
        n_clusters=optimal_n_clusters,
        linkage=config['clustering_params']['hclust_linkage']
    )
    hclust_labels = hclust.fit_predict(embedded_coords)
    hclust_silhouette = silhouette_score(embedded_coords, hclust_labels)
    print(f"Hierarchical ({config['clustering_params']['hclust_linkage']}) Silhouette Score: {hclust_silhouette:.3f}")

    results = {
        'kmeans_labels': kmeans_labels + 1,
        'hclust_labels': hclust_labels + 1,
        'optimal_n_clusters': optimal_n_clusters,
        'kmeans_silhouette': kmeans_silhouette,
        'hclust_silhouette': hclust_silhouette
    }
    
    return results

def select_cluster_representatives(embedded_coords: np.ndarray, labels: np.ndarray) -> List[int]:
    print("Selecting cluster representatives...")
    representative_frames = []
    for i in sorted(set(labels)):
        cluster_points = embedded_coords[labels == i]
        cluster_center = np.mean(cluster_points, axis=0)
        distances = np.linalg.norm(cluster_points - cluster_center, axis=1)
        norm_points = np.linalg.norm(cluster_points, axis=1)
        norm_center = np.linalg.norm(cluster_center)
        cosines = np.sum(cluster_points * cluster_center, axis=1) / (norm_points * norm_center + 1e-10)
        idx = np.where(labels == i)[0]
        min_dist_idx = idx[np.argmin(distances)]
        max_cosine_idx = idx[np.argmax(cosines)]
        representative_idx = max_cosine_idx if cosines[np.argmax(cosines)] > 0.95 else min_dist_idx  # Reverted to 0.95
        representative_frames.append(representative_idx)
        print(f"Cluster {i}: Distance-based idx={min_dist_idx}, Cosine-based idx={max_cosine_idx}, Chosen idx={representative_idx}, Max cosine={cosines[np.argmax(cosines)]:.3f}")
    
    return representative_frames

def compute_representative_rmsd(traj: md.Trajectory, binding_site_atoms: np.ndarray, labels: np.ndarray, representative_frames: List[int]) -> Dict[int, float]:
    print("Computing RMSD of representatives to their clusters...")
    binding_site_traj = traj.atom_slice(binding_site_atoms)
    rmsd_results = {}
    
    for i, rep_frame_idx in enumerate(representative_frames, 1):
        cluster_frames = np.where(labels == i)[0]
        if len(cluster_frames) > 1:
            rep_traj = binding_site_traj[rep_frame_idx]
            cluster_traj = binding_site_traj[cluster_frames]
            rmsd_values = md.rmsd(cluster_traj, rep_traj)
            avg_rmsd = np.mean(rmsd_values)
            rmsd_results[i] = avg_rmsd
        else:
            rmsd_results[i] = 0.0
        print(f"Cluster {i} representative (frame {rep_frame_idx}): Avg RMSD to cluster = {rmsd_results[i]:.3f} Å")
    
    return rmsd_results

def create_visualizations(
    embedded_coords: np.ndarray,
    kmeans_labels: np.ndarray,
    hclust_labels: np.ndarray,
    rmsd_matrix: np.ndarray,
    representative_frames: List[int],
    config: Dict
):
    print("Creating visualizations...")
    os.makedirs(config['output_dir'], exist_ok=True)
    
    plt.figure(figsize=(12, 8))
    scatter = plt.scatter(
        embedded_coords[:, 0],
        embedded_coords[:, 1],
        c=kmeans_labels,
        cmap='tab20',
        alpha=0.6
    )
    plt.scatter(
        embedded_coords[representative_frames, 0],
        embedded_coords[representative_frames, 1],
        c='red',
        marker='*',
        s=200,
        label='Centroid Frames'
    )
    plt.title('UMAP Clustering of Binding Site Conformations (K-means, 2D Projection)')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.colorbar(scatter, label='Cluster')
    plt.legend()
    plt.savefig(f"{config['output_dir']}/umap_kmeans_clusters.png", dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(12, 8))
    scatter_hclust = plt.scatter(
        embedded_coords[:, 0],
        embedded_coords[:, 1],
        c=hclust_labels,
        cmap='tab20',
        alpha=0.6
    )
    plt.scatter(
        embedded_coords[representative_frames, 0],
        embedded_coords[representative_frames, 1],
        c='red',
        marker='*',
        s=200,
        label='Centroid Frames'
    )
    plt.title(f"UMAP Clustering of Binding Site Conformations (Hierarchical, {config['clustering_params']['hclust_linkage'].capitalize()}, 2D Projection)")
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.colorbar(scatter_hclust, label='Cluster')
    plt.legend()
    plt.savefig(f"{config['output_dir']}/umap_hclust_clusters.png", dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(10, 8))
    sns.heatmap(rmsd_matrix, cmap='viridis')
    plt.title('RMSD Matrix of Binding Site')
    plt.xlabel('Frame')
    plt.ylabel('Frame')
    plt.savefig(f"{config['output_dir']}/rmsd_matrix.png", dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(10, 6))
    unique_labels = sorted(set(hclust_labels))  # Use hclust_labels for consistency
    counts = [np.sum(hclust_labels == label) for label in unique_labels]
    sns.barplot(x=unique_labels, y=counts)
    plt.title('Hierarchical (Ward) Cluster Size Distribution')
    plt.xlabel('Cluster ID')
    plt.ylabel('Number of Frames')
    plt.savefig(f"{config['output_dir']}/hclust_cluster_sizes.png", dpi=300, bbox_inches='tight')
    plt.close()

def save_results(
    traj: md.Trajectory,
    labels: np.ndarray,
    representative_frames: List[int],
    binding_site_atoms: np.ndarray,
    config: Dict,
    clustering_results: Dict
):
    print("Saving results...")
    os.makedirs(config['output_dir'], exist_ok=True)
    
    for i, frame in enumerate(representative_frames):
        traj[frame].save_pdb(f"{config['output_dir']}/cluster_{i + 1}_representative.pdb")
    
    rmsd_results = compute_representative_rmsd(traj, binding_site_atoms, labels, representative_frames)
    
    with open(f"{config['output_dir']}/analysis_summary.txt", 'w') as f:
        f.write("UMAP Analysis Summary\n")
        f.write("===================\n\n")
        f.write(f"Total frames analyzed: {len(labels)}\n")
        f.write(f"Stride used: {config['stride']}\n")
        f.write(f"Optimal number of clusters: {clustering_results['optimal_n_clusters']}\n")
        f.write(f"UMAP dimensions: {CONFIG['umap_params']['n_components']}\n")
        f.write(f"Hierarchical clustering linkage: {CONFIG['clustering_params']['hclust_linkage']}\n")
        f.write(f"K-means Silhouette Score: {clustering_results['kmeans_silhouette']:.3f}\n")
        f.write(f"Hierarchical Silhouette Score: {clustering_results['hclust_silhouette']:.3f}\n\n")
        
        f.write("Cluster statistics (Hierarchical):\n")
        for i in sorted(set(labels)):
            cluster_size = np.sum(labels == i)
            representative_frame = representative_frames[i - 1]
            avg_rmsd = rmsd_results.get(i, 0.0)
            f.write(f"Representative Cluster {i} from frame {representative_frame}\n")
            f.write(f"Cluster {i}: {cluster_size} frames ({cluster_size/len(labels)*100:.1f}%)\n")
            f.write(f"Representative frame index: {representative_frame}\n")
            f.write(f"Average RMSD to cluster: {avg_rmsd:.3f} Å\n\n")

def pca(embedded_coords: np.ndarray, labels: np.ndarray, cluster_alg: str = ''):
    # Run PCA
    pca = PCA(n_components=3)
    pca_coords = pca.fit_transform(embedded_coords) 

    # Plot PCA with hierarchical cluster labels
    plt.figure(figsize=(12, 8))
    scatter = plt.scatter(
        pca_coords[:, 0], pca_coords[:, 1],
        c=labels, 
        cmap='tab20', alpha=0.6
    )
    plt.title('PCA with Hierarchical Clustering')
    plt.xlabel('PCA 1')
    plt.ylabel('PCA 2')
    plt.colorbar(scatter, label='Cluster')
    plt.savefig(f"{CONFIG['output_dir']}/pca_sklearn_clusters_{cluster_alg}.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    try:
        traj, binding_site_atoms = load_trajectory(CONFIG)
        rmsd_matrix = calculate_rmsd_matrix(traj, binding_site_atoms)
        coords = traj.atom_slice(binding_site_atoms).xyz.reshape(traj.n_frames, -1)
        
        best_umap_params = optimize_umap_parameters(coords, CONFIG)
        print(f"Best UMAP parameters: {best_umap_params}")
        
        umap = UMAP(
            n_components=CONFIG['umap_params']['n_components'],
            **best_umap_params,
            random_state=CONFIG['random_seed']
        )
        embedded_coords = umap.fit_transform(coords)
        
        clustering_results = perform_clustering(embedded_coords, CONFIG)
        # Use Ward's labels (hclust_labels) for downstream analysis
        representative_frames = select_cluster_representatives(embedded_coords, clustering_results['hclust_labels'])
        
        create_visualizations(
            embedded_coords,
            clustering_results['kmeans_labels'],  # Still generate K-means plot for comparison
            clustering_results['hclust_labels'],
            rmsd_matrix,
            representative_frames,
            CONFIG
        )
        
        save_results(
            traj,
            clustering_results['hclust_labels'],  # Use hclust_labels
            representative_frames,
            binding_site_atoms,
            CONFIG,
            clustering_results
        )

        # PCA
        pca(coords, clustering_results['hclust_labels'], 'hclust')
        pca(coords, clustering_results['kmeans_labels'], 'kmeans')
        
        # 2D Temporal analysis plot (Cluster ID vs. Time)
        plt.figure(figsize=(12, 6))
        plt.scatter(range(len(clustering_results['hclust_labels'])), clustering_results['hclust_labels'], 
                   c=clustering_results['hclust_labels'], cmap='tab20', alpha=0.5)
        plt.title('Cluster Assignment Over Time (Hierarchical)')
        plt.xlabel('Frame Number (~0.1 ns/frame)')
        plt.ylabel('Cluster ID')
        plt.savefig(f"{CONFIG['output_dir']}/hclust_cluster_timeline.png", dpi=300)
        plt.close()

        plt.figure(figsize=(12, 6))
        plt.scatter(range(len(clustering_results['kmeans_labels'])), clustering_results['kmeans_labels'], 
                   c=clustering_results['kmeans_labels'], cmap='tab20', alpha=0.5)
        plt.title('Cluster Assignment Over Time (Kmeans)')
        plt.xlabel('Frame Number (~0.1 ns/frame)')
        plt.ylabel('Cluster ID')
        plt.savefig(f"{CONFIG['output_dir']}/kmeans_cluster_timeline.png", dpi=300)
        plt.close()
        
        # 3D Temporal analysis plot (UMAP 1-3 vs. Time)
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        scatter = ax.scatter(embedded_coords[:, 0], embedded_coords[:, 1], embedded_coords[:, 2], 
                            c=np.arange(len(embedded_coords)), cmap='viridis', alpha=0.5)
        ax.set_title('3D UMAP Timeline (Colored by Time)')
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_zlabel('UMAP 3')
        plt.colorbar(scatter, label='Frame Number (~0.1 ns/frame)')
        plt.savefig(f"{CONFIG['output_dir']}/umap_3d_timeline.png", dpi=300)
        plt.close()
        
        print("Analysis completed successfully!")

        # Save results
        labels_df = pd.DataFrame({
            'Frame': range(len(clustering_results['hclust_labels'])),
            'Cluster': clustering_results['hclust_labels']
        })
        labels_df.to_csv(f"{CONFIG['output_dir']}/hclust_labels.csv", index=False)

        labels_df = pd.DataFrame({
            'Frame': range(len(clustering_results['kmeans_labels'])),
            'Cluster': clustering_results['kmeans_labels']
        })
        labels_df.to_csv(f"{CONFIG['output_dir']}/kmeans_labels.csv", index=False)

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
