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
    nchains = len(wb_chains)

    fig, axes = plt.subplots(nchains, 1, figsize=(12, 4 * nchains), sharey=True)
    for i, chain in enumerate(wb_chains):
        chain_p1 = df_ext[df_ext['sele1_chain'] == chain].copy()
        chain_p1['sele1_resid_adj_str'] = chain_p1['sele1_resid_adj'].astype(str)
        chain_p1 = chain_p1['sele1_resid_adj_str'].value_counts().sort_index()
        axes[i].bar(chain_p1.index, chain_p1.values, color='skyblue', edgecolor='black')
        chain_p2 = df_ext[df_ext['sele2_chain'] == chain].copy()
        chain_p2['sele2_resid_adj_str'] = chain_p2['sele2_resid_adj'].astype(str)
        chain_p2 = chain_p2['sele2_resid_adj_str'].value_counts().sort_index()
        axes[i].bar(chain_p2.index, chain_p2.values, color='salmon', edgecolor='black', alpha=0.7)
        axes[i].legend([f'As atom1 in interaction', 'As atom2 in interaction'])
        axes[i].set_xlabel('Residue (adjusted)')
        axes[i].set_ylabel('Water Bridge Count')
        axes[i].set_title(f'Water Bridges Involving Chain {chain}')
        axes[i].grid(True, alpha=0.3)

    plt.suptitle(f'Water Bridges between {group1} and {group2} by Chain and Residue', fontsize=16, y=0.99)
    plt.subplots_adjust(top=0.3)
    plt.tight_layout()
    plt.savefig(outfiles['wb_per_chain'])
    plt.close()

        


if __name__ == "__main__":
    main()