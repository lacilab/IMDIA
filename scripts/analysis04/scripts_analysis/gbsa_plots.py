#!/usr/bin/env python3

from io import StringIO
from pathlib import Path
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from scripts.utils.helpers import setup_logger

def get_delta_total (logger, energyfile: Path):
    with open(energyfile, 'r') as f:
        lines = f.readlines()

    # Find the index of the line containing "Delta Energy Terms"
    start_idx = None
    for i, line in enumerate(lines):
        if "Delta Energy Terms" in line:
            start_idx = i + 1  # Start from the next line
            break

    # If "Delta Energy Terms" not found
    if start_idx is None:
        logger.error("Section 'Delta Energy Terms' not found in the file. The GBSA results may be incomplete. Try rerunning --analysis --rules gbsa.")
        return
    
    # Extract the lines after "Delta Energy Terms"
    data_lines = lines[start_idx:]
    data_string = ''.join(data_lines)
    df = pd.read_csv(StringIO(data_string))
    return df

def load_delta_decomp(logger, energyfile: Path):
    with open(energyfile, 'r') as f:
            lines = f.readlines()
    # Find the start and end indices for the section
    start_idx = None
    end_idx = None
    for i, line in enumerate(lines):
        line_clean = line.strip()
        if "DELTAS:" in line_clean:
            start_idx = i + 2  # Skip "DELTAS:" and "Total Decomposition Contribution (TDC)"
            continue
        if start_idx is not None and "Sidechain Decomposition Contribution" in line_clean:
            end_idx = i  
            break

    # If the section is not found
    if start_idx is None or end_idx is None:
        logger.error(f"Section 'DELTAS:' not found in {energyfile}. The GBSA results may be incomplete. Try rerunning --analysis --rules gbsa.")
        return

    # Extract the lines between start_idx and end_idx
    data_lines = lines[start_idx:end_idx]
    data_string = ''.join(data_lines)

    df = pd.read_csv(StringIO(data_string))
    return df

def main(logfile: str, energy_total: Path, energy_decomp: Path, offsets: dict, outdir: Path):
    """ 
    Plot GBSA result.

    Arguments:
        logfile        : Path to the log file.
        energy_total   : Path to .csv file with total energy MM/GBSA results, e.g. rep1/analysis/gbsa/FINAL_ENERGY.csv
        energy_decomp  : Path to .csv file with MM/GBSA results decomposed by residue, e.g. rep1/analysis/gbsa/FINAL_ENERGY_DECOMP.csv
        offsets        : Dictionary specifying residue numbering offsets per chain.
        outdir         : Output directory for MM/GBSA results (e.g. rep1/analysis/gbsa)
    """
    # Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug("Plotting GBSA results")

    # Expected output files
    filenames = [
        'gbsa_full_totalenergy.png',
        'gbsa_res_totalenergy.png'

    ]
    outdir = Path(outdir)
    outfiles = [outdir / n for n in filenames]

    # Check files exist and are .csv
    if not energy_total.is_file() or energy_total.suffix != '.csv':
        logger.error(f"Energy total file {energy_total} does not exist or is not a .csv file.")
        return
    if not energy_decomp.is_file() or energy_decomp.suffix != '.csv':
        logger.error(f"Energy decomposition file {energy_decomp} does not exist or is not a .csv file.")
        return
    
    ##### TOTAL ENERGY ANALYSIS #####
    
    df_total = get_delta_total(logger, energy_total)
    
    # Get stats for total energy
    components = ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL', 'EGB', 'ESURF', 'GGAS', 'GSOLV', 'TOTAL']
    df_total_desc = df_total[components].describe().filter(items = ['mean', 'std'], axis=0).transpose()
    means = df_total_desc['mean'].values
    stds = df_total_desc['std'].values

    # Define groups based on GBSA relationships
    group_map = {}
    for c in components:
        if c in ['BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL']:
            group_map[c] = 'GGAS'
        elif c in ['EGB', 'ESURF']:
            group_map[c] = 'GSOLV'
        elif c == 'TOTAL':
            group_map[c] = 'TOTAL'
        else:  # GGAS, GSOLV if present
            group_map[c] = c

    colors = {'GGAS': 'skyblue', 'GSOLV': 'lightgreen', 'TOTAL': 'coral'}
    bar_colors = [colors[group_map[c]] for c in components]

    # Plot total energy
    plt.figure(figsize=(8, 8))
    plt.bar(components, means, yerr=stds, capsize=5, color=bar_colors, edgecolor='black')
    plt.ylabel('Energy (kcal/mol)')
    plt.xlabel('Energy Component')
    plt.xticks(rotation=90)
    plt.title("GBSA Total Energy")
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(outdir / "gbsa_full_totalenergy.png")
    plt.close()


    ##### ENERGY DECOMPOSITION BY RESIDUE ANALYSIS #####
    df_decomp = load_delta_decomp(logger, energy_decomp)

    # Map residue ids
    split_cols = df_decomp["Residue"].str.split(":", expand=True)
    df_decomp["compound"] = split_cols[0]
    df_decomp["chain"] = split_cols[1]
    df_decomp["resnm"] = split_cols[2]
    df_decomp["resid"] = split_cols[3]
    df_decomp["resid_adj"] = (df_decomp['resid'].astype(int) + df_decomp['chain'].map(
        lambda chain: offsets.get(chain, {}).get('offset', 0)
    ).astype(int)) # Adjust resid by chain offset and convert back to string
    # Adjust chain id, new->old
    df_decomp["chain"] = df_decomp["chain"].map(lambda x: offsets.get(x, {}).get("id_old", x))

    # Get stats by residue for total energy
    # Group by residue identifiers
    grouped = df_decomp.groupby(['Residue', 'compound', 'chain', 'resnm', 'resid', 'resid_adj'])
    
    # Compute mean, std, and sem for energy columns
    df_decomp_desc = pd.concat([
        grouped['TOTAL'].mean().rename("TOTAL_mean"),
        grouped['TOTAL'].std().rename('TOTAL_std'),
    ], axis=1)
    
    # Reset index to make Residue, compound, chain, resnm, resid regular columns
    df_decomp_desc = df_decomp_desc.reset_index()

    # Plot total energy by residue. One plot per chain, in alphabetical order
    named_chains = sorted([c for c in df_decomp_desc['chain'].unique() if len(c.strip()) > 0])
    nchains = len(named_chains)

    # Create GridSpec: nchains rows, 2 columns - plots use col 0 only (narrow for long chains)
    fig, axes = plt.subplots(nchains, 1, figsize=(12, 4*nchains), sharey=True)

    for i, chain in enumerate(named_chains):
        # Sort by residue and convert to string
        df_chain_sorted = df_decomp_desc[df_decomp_desc['chain'] == chain].sort_values('resid_adj')
        x = df_chain_sorted['resid_adj'].astype(str)  # categorical labels
        y = df_chain_sorted['TOTAL_mean']
        yerr = df_chain_sorted['TOTAL_std']

        axes[i].bar(x, y, yerr=yerr, capsize=5, color='skyblue', edgecolor='black')
        
        axes[i].axhline(0, color='black', linewidth=0.5)
        axes[i].set_xticks(x)
        axes[i].set_xticklabels(df_chain_sorted['resid_adj'], rotation=90)
        axes[i].set_xlabel('Residue (adj)')
        axes[i].set_ylabel('Energy (kcal/mol)')
        axes[i].grid(True, axis='y', linestyle='--', alpha=0.7)
        axes[i].set_title(f"Chain {chain}")

    plt.suptitle("GBSA Total Energy per Residue", fontsize=16, y=0.99)
    plt.subplots_adjust(top=0.3)
    plt.tight_layout()
    plt.savefig(outdir / "gbsa_res_totalenergy.png")
    plt.close()

    # Check all expected outputs were generated
    missing = [f.name for f in outfiles if not f.is_file()]
    if missing:
        logger.warning("GBSA plotting did not complete successfully. Some expected outputs are missing.")
        logger.warning(f"Missing {len(missing)} file(s): {', '.join(missing)}")
        logger.warning(f"Check log {logfile}")
    else:
        logger.info(f"GBSA plotting completed successfully. Output files are in {outdir}.")


if __name__ == "__main__":
    main()