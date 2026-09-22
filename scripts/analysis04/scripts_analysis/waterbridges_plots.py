#!/usr/bin/env python3

import pandas as pd
import numpy as np
import MDAnalysis as mda
from matplotlib import pyplot as plt
from scripts.utils.helpers import setup_logger

def main(logfile, wb_csv, wb_ext_csv, outfiles, group1, group2):
    """
    Plot water-bridge analysis results between groups.

    Arguments:
        logfile      : Path to the log file.
        wb_csv       : Path to .csv file containing raw water-bridge results.
        wb_ext_csv   : Path to .csv file with extended water bridge results.
        outfiles     : Dictionary of the files to be generated, assigned to labels/keys.
        group1       : First atom group for hydrogen bond analysis.
        group2       : Second atom group for hydrogen bond analysis.
    """

    logger = setup_logger("analysis", logfile)
    logger.debug("Plotting water bridge analysis results.")

    # --- Load data safely ---
    try:
        df = pd.read_csv(wb_csv)
    except Exception:
        logger.warning(f"Failed to read {wb_csv}, using empty dataframe.")
        df = pd.DataFrame()

    try:
        df_ext = pd.read_csv(wb_ext_csv)
    except Exception:
        logger.warning(f"Failed to read {wb_ext_csv}, using empty dataframe.")
        df_ext = pd.DataFrame()

    # ============================================================
    # 1. WATER BRIDGES PER FRAME
    # ============================================================

    plt.figure(figsize=(10, 6))

    if not df.empty and 'time' in df.columns:
        time_counts = df['time'].value_counts().sort_index()
        plt.plot(time_counts.index, time_counts.values)
    else:
        logger.debug("No data available for WB per frame plot.")
        plt.text(0.5, 0.5, "No water bridges detected",
                 ha='center', va='center', transform=plt.gca().transAxes)

    plt.xlabel('Time (ps)')
    plt.ylabel('Water Bridge Count')
    plt.title(f'Water Bridge Count between {group1} and {group2}')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outfiles['wb_per_frame'])
    plt.close()

    # ============================================================
    # 2. WATER BRIDGES PER CHAIN / RESIDUE
    # ============================================================

    required_cols = {
        'sele1_chain', 'sele2_chain',
        'sele1_resid_adj', 'sele2_resid_adj'
    }

    all_p1 = pd.Series(dtype=int)
    all_p2 = pd.Series(dtype=int)

    if not df_ext.empty and required_cols.issubset(df_ext.columns):

        sele1_chains = [c for c in df_ext['sele1_chain'].unique() if c != 'UNK']
        sele2_chains = [c for c in df_ext['sele2_chain'].unique() if c != 'UNK']
        wb_chains = sorted(set(sele1_chains + sele2_chains))

        p1_list = []
        p2_list = []

        for chain in wb_chains:

            # Atom1 contributions
            chain_p1 = df_ext[df_ext['sele1_chain'] == chain]
            if not chain_p1.empty:
                labels = f"{chain}:" + chain_p1['sele1_resid_adj'].astype(str)
                p1_list.append(labels.value_counts())

            # Atom2 contributions
            chain_p2 = df_ext[df_ext['sele2_chain'] == chain]
            if not chain_p2.empty:
                labels = f"{chain}:" + chain_p2['sele2_resid_adj'].astype(str)
                p2_list.append(labels.value_counts())

        if p1_list:
            all_p1 = pd.concat(p1_list)

        if p2_list:
            all_p2 = pd.concat(p2_list)

    else:
        logger.debug("No valid data for WB per chain plot.")

    # --- Align indices safely ---
    shared_labels = all_p1.index.union(all_p2.index)
    all_p1 = all_p1.reindex(shared_labels, fill_value=0)
    all_p2 = all_p2.reindex(shared_labels, fill_value=0)

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(max(4 * len(shared_labels)//5, 6), 6))

    if len(shared_labels) > 0:
        ax.bar(shared_labels, all_p1.values,
               color='skyblue', edgecolor='black',
               label='As atom1 in interaction')

        ax.bar(shared_labels, all_p2.values,
               color='salmon', edgecolor='black',
               alpha=0.7,
               label='As atom2 in interaction')
    else:
        ax.text(0.5, 0.5, "No water bridges detected",
                ha='center', va='center', transform=ax.transAxes)

    ax.set_xlabel('Chain:Residue (adjusted)')
    ax.set_ylabel('Water Bridge Count')
    ax.set_title('Water Bridge Count')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(outfiles['wb_per_chain'])
    plt.close()

    logger.debug("Finished plotting water bridge analysis.")