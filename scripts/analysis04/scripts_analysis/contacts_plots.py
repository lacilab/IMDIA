#!/usr/bin/env python3

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
    df_pivot = df_pairs_ext.pivot_table(index='sele1_label', columns='sele2_label', values='ratio')
    df_pivot = df_pivot.reindex(
        index=sorted(df_pivot.index, key=lambda x: label_sort_key(str(x))),
        columns=sorted(df_pivot.columns, key=lambda x: label_sort_key(str(x)))
    )
    # Dynamic size based on number of residues
    nx = len(df_pivot.columns)
    ny = len(df_pivot.index)
    fig = plt.figure(figsize=(max(10, nx * 0.25), max(8, ny * 0.25)))
    sns.heatmap(df_pivot, cmap='rocket_r', linewidths=0.5, 
                linecolor='lightgray', cbar_kws={'label': 'Count / Number of frames', "shrink": 0.4})
    plt.title(f'Contact Count (< 5 Å) between {group1} and {group2} by Residue Pair', fontsize=16, y=1.02)
    plt.xticks(rotation=-90)
    plt.ylabel(f'{group1} Residue (adj)', fontsize=14)
    plt.xlabel(f'{group2} Residue (adj)', fontsize=14)
    plt.tight_layout()
    plt.savefig(outfiles['contacts_heatmap_count'])
    plt.close()

    ##### Heatmap - Average distance per residue pair ######
    df_pivot = df_pairs_ext.pivot_table(index='sele1_label', columns='sele2_label', values='distance_avg')
    df_pivot = df_pivot.reindex(
        index=sorted(df_pivot.index, key=lambda x: label_sort_key(str(x))),
        columns=sorted(df_pivot.columns, key=lambda x: label_sort_key(str(x)))
    )
    # Dynamic size based on number of residues
    nx = len(df_pivot.columns)
    ny = len(df_pivot.index)
    fig = plt.figure(figsize=(max(10, nx * 0.25), max(8, ny * 0.25)))
    sns.heatmap(df_pivot, cmap='rocket', linewidths=0.5, 
                linecolor='lightgray', cbar_kws={'label': 'Average Distance (Å)', "shrink": 0.4})
    plt.title(f'Average Contact Distance between {group1} and {group2} by Residue Pair', fontsize=16, y=1.02)
    plt.xticks(rotation=-90)
    plt.ylabel(f'{group1} Residue (adj)', fontsize=14)
    plt.xlabel(f'{group2} Residue (adj)', fontsize=14)
    plt.tight_layout()
    plt.savefig(outfiles['contacts_heatmap_distance'])
    plt.close()

if __name__ == "__main__":
    main()