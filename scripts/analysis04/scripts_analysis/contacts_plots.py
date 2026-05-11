#!/usr/bin/env python3

import math
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from scripts.utils.helpers import setup_logger


def label_sort_key(label):
    if ':' in label:
        chain, res = label.split(':', 1)
        try:
            return (chain, int(res))
        except:
            return (chain, res)
    else:
        try:
            return ('@', int(label))  # numbers first
        except:
            return (label, 0)         # names last

def main(logfile: str, universe: str, pairs_csv: Path, pairs_ext_csv: Path, group1: str, group2: str, outfiles: dict):
    # Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug("Parsing contact analysis results.")

    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[-1] 
    nframes = u.trajectory.frame

    df_pairs = pd.read_csv(pairs_csv, sep=',')
    df_pairs_ext = pd.read_csv(pairs_ext_csv, sep=',')


    ##### Number of contacts vs distance (sanity check) #####
    weights = np.ones(len(df_pairs)) / nframes
    plt.hist(df_pairs['distance'], bins=50, alpha=0.5, weights=weights)
    plt.xlabel('Distance (Å)')
    plt.ylabel('Count / Number of frames')
    plt.title(f'Distance Distribution of Contacts between {group1} and {group2}')
    plt.grid(True, alpha=0.3)
    plt.savefig(outfiles['contacts_distance'])
    plt.close()


    ##### Number of contacts vs time ######
    plt.hist(df_pairs['frame'], bins=100, alpha=0.5)
    plt.xlabel('Frame')
    plt.ylabel('Count')
    plt.title(f'Time Distribution of Contacts between {group1} and {group2}')
    plt.grid(True, alpha=0.3)
    plt.savefig(outfiles['contacts_time'])
    plt.close()


    ##### Heatmap - Number of contacts per residue pair ######
    df_pairs_ext['sele1_label'] = (df_pairs_ext['sele1_chain'].astype(str) + ":" + 
                              df_pairs_ext['sele1_resid_adj'].astype(str))
    df_pairs_ext['sele2_label'] = (df_pairs_ext['sele2_chain'].astype(str) + ":" + 
                              df_pairs_ext['sele2_resid_adj'].astype(str))
    df_pivot_ratio = df_pairs_ext.pivot_table(index='sele1_label', columns='sele2_label', values='ratio')
    df_pivot_ratio = df_pivot_ratio.reindex(
        index=sorted(df_pivot_ratio.index, key=lambda x: label_sort_key(str(x))),
        columns=sorted(df_pivot_ratio.columns, key=lambda x: label_sort_key(str(x)))
    )
    # Dynamic size based on number of residues
    nx = len(df_pivot_ratio.columns)
    ny = len(df_pivot_ratio.index)
    fig = plt.figure(figsize=(max(10, nx * 0.25), max(8, ny * 0.25)))
    sns.heatmap(df_pivot_ratio, cmap='rocket_r', linewidths=0.5, 
                linecolor='lightgray', cbar_kws={'label': 'Count / Number of frames', "shrink": 0.4})
    plt.title(f'Contact Count (< 5 Å) between {group1} and {group2} by Residue Pair', fontsize=16, y=1.02)
    plt.xticks(rotation=-90)
    plt.ylabel(f'{group1} Residue (adj)', fontsize=14)
    plt.xlabel(f'{group2} Residue (adj)', fontsize=14)
    plt.tight_layout()
    plt.savefig(outfiles['contacts_heatmap_count'])
    plt.close()

    ##### Heatmap - Average distance per residue pair ######
    df_pivot_dist = df_pairs_ext.pivot_table(index='sele1_label', columns='sele2_label', values='distance_avg')
    df_pivot_dist = df_pivot_dist.reindex(
        index=sorted(df_pivot_dist.index, key=lambda x: label_sort_key(str(x))),
        columns=sorted(df_pivot_dist.columns, key=lambda x: label_sort_key(str(x)))
    )
    # Dynamic size based on number of residues
    nx = len(df_pivot_dist.columns)
    ny = len(df_pivot_dist.index)
    fig = plt.figure(figsize=(max(10, nx * 0.25), max(8, ny * 0.25)))
    sns.heatmap(df_pivot_dist, cmap='rocket', linewidths=0.5, 
                linecolor='lightgray', cbar_kws={'label': 'Mean Distance (Å)', "shrink": 0.4})
    plt.title(f'Mean Contact Distance between {group1} and {group2} by Residue Pair', fontsize=16, y=1.02)
    plt.xticks(rotation=-90)
    plt.ylabel(f'{group1} Residue (adj)', fontsize=14)
    plt.xlabel(f'{group2} Residue (adj)', fontsize=14)
    plt.tight_layout()
    plt.savefig(outfiles['contacts_heatmap_distance'])
    plt.close()

    ##### COMBINED BUBBLE PLOT - Count and Average Distance ######

    # Grid coordinates
    x_labels = df_pivot_ratio.columns.tolist()
    y_labels = df_pivot_ratio.index.tolist()
    x = np.arange(len(x_labels))
    y = np.arange(len(y_labels))
    xx, yy = np.meshgrid(x, y)

    # Flatten data
    ratio_vals = df_pivot_ratio.to_numpy().ravel()
    dist_vals = df_pivot_dist.to_numpy().ravel()

    min_ratio = math.floor(np.nanmin(ratio_vals))
    max_ratio = math.ceil(np.nanmax(ratio_vals))
    step = math.ceil((max_ratio - min_ratio) / 4)

    # Size scaling
    size_scale = 1400 / max_ratio
    sizes = ratio_vals * size_scale
    padding = 0.6

    fig, ax = plt.subplots(figsize=(max(10, len(x_labels) * 0.5), max(8, len(y_labels) * 0.5)))

    sc = ax.scatter(xx, yy, s=sizes, c=dist_vals, cmap='rocket', alpha=0.9,
                    edgecolors='lightgray', linewidths=0.5)

    # EXPLICIT tight limits (removes big margins)
    ax.set_xlim(-0.5 - padding, len(x_labels)-0.5 + padding)
    ax.set_ylim(-0.5 - padding, len(y_labels)-0.5 + padding)

    # Set ticks FIRST
    ax.set_xticks(np.arange(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=-90)
    ax.set_yticks(np.arange(len(y_labels)))
    ax.set_yticklabels(y_labels)

    # Invert Y-axis
    ax.invert_yaxis()


    # Labels/title
    ax.set_xlabel(f'{group2} Residue (adj)', fontsize=14)
    ax.set_ylabel(f'{group1} Residue (adj)', fontsize=14)
    ax.set_title(
        f'Contacts between {group1} and {group2}',
        fontsize=20, y=1.02
    )

    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, shrink=0.4)
    cbar.set_label('Mean Distance (Å)', fontsize=14)


    # Size legend (adjust bbox_to_anchor if needed)
    size_handles = []
    for r in np.arange(min_ratio + step, max_ratio + step, step):
        size_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                    markerfacecolor='lightgray',
                                    markersize=np.sqrt(r * size_scale),
                                    markeredgecolor='black', label=f'{r:.1f}'))
    fig.legend(handles=size_handles, title='Ratio (Count / Number of frames)', 
            loc='center left', bbox_to_anchor=(0.82, 0.8),
            ncol=1, frameon=True, fontsize=12,
            labelspacing=2, handletextpad=1.2, borderpad=1.5)

    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outfiles['contacts_bubble'], bbox_inches='tight', dpi=300, pad_inches=0.1)
    plt.close()

    ##### COMBINED BUBBLE PLOT - FILTERED - Count and Average Distance ######
    # Filter to pairs with ratio in the top 75%
    threshold = round(np.nanpercentile(ratio_vals, 75),0)
    filter_mask = ratio_vals >= threshold
    x_filtered = xx.ravel()[filter_mask]
    y_filtered = yy.ravel()[filter_mask]
    sizes_filtered = sizes[filter_mask]
    dist_filtered = dist_vals[filter_mask]


    # Full grid ticks/labels (no empty rows)
    unique_x = np.sort(np.unique(x_filtered).astype(int))
    unique_y = np.sort(np.unique(y_filtered).astype(int))

    # Remap to dense coordinates [0,1,2,...] instead of sparse [2,5,7,12,...]
    x_remap = {old: new for new, old in enumerate(unique_x)}
    y_remap = {old: new for new, old in enumerate(unique_y)}

    # Transform filtered coordinates
    x_dense = np.array([x_remap[int(xi)] for xi in x_filtered])
    y_dense = np.array([y_remap[int(yi)] for yi in y_filtered])

    fig, ax = plt.subplots(figsize=(max(10, len(x_dense) * 0.3), max(8, len(y_dense) * 0.3)))
    sc = ax.scatter(x_dense, y_dense, s=sizes_filtered, c=dist_filtered,
                    cmap='rocket', alpha=0.9, edgecolors='lightgray', linewidths=0.5)

    # 1. MANUAL padding (margins conflicts with set_xlim/ylim)
    ax.set_xlim(-0.5-padding, len(unique_x)-0.5+padding)
    ax.set_ylim(-0.5-padding, len(unique_y)-0.5+padding)

    # 2. Set ticks FIRST (before inversion)
    ax.set_xticks(np.arange(len(unique_x)))
    ax.set_xticklabels([x_labels[i] for i in unique_x], rotation=-90)
    ax.set_yticks(np.arange(len(unique_y)))
    ax.set_yticklabels([y_labels[i] for i in unique_y])

    # 3. Invert Y-axis AFTER ticks
    ax.invert_yaxis()

    # 4. Labels/title
    ax.set_xlabel(f'{group2} Residue (adj)', fontsize=14)
    ax.set_ylabel(f'{group1} Residue (adj)', fontsize=14)
    ax.set_title(
        f'Contacts between {group1} and {group2}\n'
        f'Filtered to pairs with ratio >= {threshold:.2f} (top 75%)',
        fontsize=20, y=1.02
    )

    # 5. Colorbar
    cbar = plt.colorbar(sc, ax=ax, shrink=0.4)
    cbar.set_label('Mean Distance (Å)', fontsize=14)

    # 6. Size legend (moved LEFT side to avoid colorbar overlap)
    step = math.ceil((max_ratio - threshold) / 3)
    size_handles = []
    for r in np.arange(threshold, max_ratio + step, step):
        size_handles.append(
            plt.Line2D([0], [0], marker='o', color='w', 
                    markerfacecolor='lightgray', 
                    markersize=np.sqrt(r * size_scale),
                    markeredgecolor='black', 
                    label=f'{r:.1f}')
        )

    fig.legend(handles=size_handles, title='Ratio (Count / Number of frames)', 
            loc='center left', bbox_to_anchor=(0.82, 0.8),
            ncol=1, frameon=True, fontsize=12,
            labelspacing=2, handletextpad=1.2, borderpad=1.5)

    # 7. Layout LAST (no conflicts)
    plt.tight_layout()
    plt.grid(True, alpha=0.3)
    plt.savefig(outfiles['contacts_bubble_filtered'], bbox_inches='tight', dpi=300)
    plt.close()


if __name__ == "__main__":
    main()