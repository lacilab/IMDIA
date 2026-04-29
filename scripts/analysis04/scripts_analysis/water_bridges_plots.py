#!/usr/bin/env python3

import pandas as pd
import numpy as np
import MDAnalysis as mda
from matplotlib import pyplot as plt
from scripts.utils.helpers import setup_logger

def main(logfile, wb_csv, wb_ext_csv, outfiles, group1, group2, chains):

    # Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug("Plotting water bridge analysis results.")

    df = pd.read_csv(wb_csv, sep=',')
    df_ext = pd.read_csv(wb_ext_csv, sep=',')

    ##### WB PER FRAME #####
    time_counts = df['time'].value_counts().sort_index()

    plt.figure(figsize=(10, 6))
    plt.plot(time_counts.index, time_counts.values, color='blue')
    plt.xlabel('Time (ps)')
    plt.ylabel('Water Bridge Count')
    plt.title(f'Water Bridge Count between {group1} and {group2}')
    plt.grid(True, alpha=0.3)
    plt.savefig(outfiles['wb_per_frame'])
    plt.close()

    ##### WB PER CHAIN AND RESIDUE #####
    # Get chain ID of chains involved in water bridges
    sele1_chains = [c for c in df_ext['sele1_chain'].unique() if c != 'UNK']
    sele2_chains = [c for c in df_ext['sele2_chain'].unique() if c != 'UNK']
    wb_chains = sorted(set(sele1_chains + sele2_chains))

    all_p1 = []
    all_p2 = []

    for chain in wb_chains:
        # Atom1 (sele1)
        chain_p1 = df_ext[df_ext['sele1_chain'] == chain].copy()
        chain_p1['label'] = f"{chain}:" + chain_p1['sele1_resid_adj'].astype(str)
        p1_counts = chain_p1['label'].value_counts().sort_index()
        all_p1.append(p1_counts)

        # Atom2 (sele2)
        chain_p2 = df_ext[df_ext['sele2_chain'] == chain].copy()
        chain_p2['label'] = f"{chain}:" + chain_p2['sele2_resid_adj'].astype(str)
        p2_counts = chain_p2['label'].value_counts().sort_index()
        all_p2.append(p2_counts)

    # Combine all chains
    all_p1 = pd.concat(all_p1, axis=0).sort_index()
    all_p2 = pd.concat(all_p2, axis=0).sort_index()

    # Re‑index both to the same set of labels
    shared_labels = all_p1.index.union(all_p2.index)
    all_p1 = all_p1.reindex(shared_labels, fill_value=0)
    all_p2 = all_p2.reindex(shared_labels, fill_value=0)

    # Plot
    fig, ax = plt.subplots(figsize=(max(4 * len(shared_labels)//5, 6), 6))
    ax.bar(
        shared_labels,
        all_p1.values,
        color='skyblue',
        edgecolor='black',
        label='As atom1 in interaction'
    )
    ax.bar(
        shared_labels,
        all_p2.values,
        color='salmon',
        edgecolor='black',
        alpha=0.7,
        label='As atom2 in interaction'
    )
    ax.set_xlabel('Chain:Residue (adjusted)')
    ax.set_ylabel('Water Bridge Count')
    ax.set_title('Water Bridge Count')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()
    plt.savefig(outfiles['wb_per_chain'])
    plt.close()

        


if __name__ == "__main__":
    main()