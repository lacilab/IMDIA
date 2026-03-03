#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
from scripts.utils.helpers import setup_logger

def get_chaingroup(chain, chaingroups):
    for group, chains in chaingroups.items():
        if chain in chains:
            return group
    return None

def main(logfile : str =None, hbonds_ext : str =None, domains : dict = {}, chaingroups: dict = {}, 
         group1: str = None, group2 : str = None, outfiles : list = []):
    # Arguments
    # parser = argparse.ArgumentParser(description='Generate plots from extended hydrogen bonds data.')
    # parser.add_argument("-i", "--incsv", type=str, help='Path to the extended hydrogen bonds CSV file')
    # parser.add_argument("-m", "--method", type=str, help='Method used for hydrogen bond calculation (e.g., "gmx, mda").')
    # parser.add_argument("-d", "--outdir", type=str, help="Output directory")

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    hbonds_df = pd.read_csv(hbonds_ext)

    ##########  Plots per group
    markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h', 'H', 'X']
    # 2 groups
    fig, axes = plt.subplots(2, 1, figsize=(10, 5*2))
    for i, group in enumerate([group1, group2]):
        # If group is a chain ID, use regex
        if len(group) == 1 and group.isalpha():
            donor_df = hbonds_df[hbonds_df['donor_chain'] == group]
            acceptor_df = hbonds_df[hbonds_df['acceptor_chain'] == group]
            group_domains = get_chaingroup(group, chaingroups)
        # If group is a predefined group (e.g., "Protein", "DNA")
        elif group == "Protein" or group == "DNA":
            donor_df = hbonds_df[hbonds_df['donor_group'] == group]
            acceptor_df = hbonds_df[hbonds_df['acceptor_group'] == group]
            group_domains = group if group in domains.keys() else None
        chains = set(donor_df['donor_chain'].unique().tolist() + acceptor_df['acceptor_chain'].unique().tolist())
        for j, chain in enumerate(chains):
            donor_chain_df = donor_df[donor_df['donor_chain'] == chain]
            acceptor_chain_df = acceptor_df[acceptor_df['acceptor_chain'] == chain]
            axes[i].scatter(donor_chain_df['donor_res_num_adjusted'], donor_chain_df['frame%'], 
                    color='red', label=f'Donor Chain {chain}', alpha=0.6, marker=markers[j % len(markers)])
            axes[i].scatter(acceptor_chain_df['acceptor_res_num_adjusted'], acceptor_chain_df['frame%'], 
                    color='blue', label=f'Acceptor Chain {chain}', alpha=0.6, marker=markers[j % len(markers)])
        axes[i].set_xlabel('Residue Number (Adjusted)')
        axes[i].set_ylabel('Frame %')
        axes[i].set_title(f'Hydrogen Bonds for {group}')
        axes[i].grid(True, alpha=0.3)
        axes[i].legend()
        # Color by domain if available
        if group_domains and group_domains in domains.keys():
            for domain in domains[group_domains]:
                axes[i].axvspan(domain['start'], domain['end'], zorder=0, alpha=0.2, color=domain['color'], label=domain['label'])

    plt.tight_layout()
    plt.savefig(outfiles["hbonds_per_group"])
    plt.close()
    logger.info(f"Plot of hydrogen bonds per group saved as {outfiles['hbonds_per_group']}")


    ########## Hbonds per chain (if groups have multiple chains)
    for group in [group1, group2]:
        if group == "Protein":
            chains = hbonds_df[hbonds_df['donor_group'] == group]['donor_chain'].unique()
            nchains = len(chains)
            fig, axes = plt.subplots(nchains, 1, figsize=(10, 5*nchains))
            if nchains == 1:
                axes = [axes]
            for i, chain in enumerate(chains):
                donor_df = hbonds_df[hbonds_df['donor_chain'] == chain]
                acceptor_df = hbonds_df[hbonds_df['acceptor_chain'] == chain]
                axes[i].scatter(donor_df['donor_res_num_adjusted'], donor_df['frame%'], 
                        color='red', label=f'Donor', alpha=0.6)
                axes[i].scatter(acceptor_df['acceptor_res_num_adjusted'], acceptor_df['frame%'], 
                        color='blue', label=f'Acceptor', alpha=0.6)
                axes[i].set_xlabel('Residue Number (Adjusted)')
                axes[i].set_ylabel('Frame %')
                axes[i].set_title(f'Hydrogen Bonds for Chain {chain}')
                axes[i].grid(True, alpha=0.3)
                axes[i].legend()
                # Color by domain if available
                chain_group = get_chaingroup(chain, chaingroups)
                if chain_group and chain_group in domains.keys():
                        for domain in domains[chain_group]:
                                axes[i].axvspan(domain['start'], domain['end'], zorder=0, alpha=0.2, color=domain['color'], label=domain['label'])    
            plt.tight_layout()
            plt.savefig(outfiles[f"hbonds_prot_chains"])
            plt.close()
            logger.info(f"Plot of protein hydrogen bonds per chain saved as {outfiles[f'hbonds_prot_chains']}")
        elif group == "DNA": 
            chains = hbonds_df[hbonds_df['donor_group'] == group]['donor_chain'].unique()
            nchains = len(chains)
            fig, axes = plt.subplots(nchains, 1, figsize=(10, 5*nchains))
            if nchains == 1:
                axes = [axes]
            for i, chain in enumerate(chains):
                donor_df = hbonds_df[hbonds_df['donor_chain'] == chain]
                acceptor_df = hbonds_df[hbonds_df['acceptor_chain'] == chain]
                sections = set(donor_df['donor_section'].unique().tolist() + acceptor_df['acceptor_section'].unique().tolist())
                for j, section in enumerate(sections):
                    donor_section_df = donor_df[donor_df['donor_section'] == section]
                    acceptor_section_df = acceptor_df[acceptor_df['acceptor_section'] == section]
                    axes[i].scatter(donor_section_df['donor_res_num_adjusted'], donor_section_df['frame%'], 
                            color='red', marker=markers[j % len(markers)],  label=f'Donor, {section}', alpha=0.6)
                    axes[i].scatter(acceptor_section_df['acceptor_res_num_adjusted'], acceptor_section_df['frame%'], 
                            color='blue', label=f'Acceptor, {section}', alpha=0.6)
                axes[i].set_xlabel('Residue Number (Adjusted)')
                axes[i].set_ylabel('Frame %')
                axes[i].set_title(f'Hydrogen Bonds for Chain {chain}')
                axes[i].grid(True, alpha=0.3)
                axes[i].legend()
                # Color by domain if available
                chain_group = get_chaingroup(chain, chaingroups)
                if chain_group and chain_group in domains.keys():
                        for domain in domains[chain_group]:
                                axes[i].axvspan(domain['start'], domain['end'], zorder=0, alpha=0.2, color=domain['color'], label=domain['label'])    
            plt.tight_layout()
            plt.savefig(outfiles[f"hbonds_dna_chains"])
            plt.close()
            logger.info(f"Plot of DNA hydrogen bonds per chain saved as {outfiles[f'hbonds_dna_chains']}")

      
if __name__ == "__main__":
    main()