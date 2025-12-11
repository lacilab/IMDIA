#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
from scripts.utils.helpers import setup_logger


def main(logfile : str =None, hbonds_ext : str =None, method : str =None, dna: bool=False, domains : list = [], outfiles : list = []):
    # Arguments
    # parser = argparse.ArgumentParser(description='Generate plots from extended hydrogen bonds data.')
    # parser.add_argument("-i", "--incsv", type=str, help='Path to the extended hydrogen bonds CSV file')
    # parser.add_argument("-m", "--method", type=str, help='Method used for hydrogen bond calculation (e.g., "gmx, mda").')
    # parser.add_argument("-d", "--outdir", type=str, help="Output directory")
# 
    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit(0)
# 
    # args = parser.parse_args()
    # hbonds_ext = args.incsv
    # outdir = args.outdir
    # method = args.method

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    hbonds_df = pd.read_csv(hbonds_ext)

    ##########  Protein acceptors and donors in the whole dimer
    donor_df = hbonds_df[hbonds_df['donor_type']=='protein']
    acceptor_df = hbonds_df[hbonds_df['acceptor_type']=='protein']

    plt.figure(figsize=(10, 6))

    # Plot donor residues in red
    plt.scatter(donor_df['donor_res_num_adjusted'], donor_df['frame%'], 
            color='red', label='Donor', alpha=0.6)

    # Plot acceptor residues in blue (if any have 3-letter codes)
    plt.scatter(acceptor_df['acceptor_res_num_adjusted'], acceptor_df['frame%'], 
            color='blue', label='Acceptor', alpha=0.6)

    plt.xlabel('Residue Number (Adjusted)')
    plt.ylabel('Frame %')
    plt.title('STAT1 Donor and Acceptor Residues')
    plt.grid(True, alpha=0.3)
    if len(domains) > 0:
        for domain in domains:
            start = domain['start']
            end = domain['end']
            plt.axvspan(start, end, zorder=0, alpha=0.2, color=domain['color'], label=domain['label'])
    plt.legend()
    plt.savefig(outfiles["hbonds_prot_domains"])
    plt.close()
    logger.info(f"Plot of protein hydrogen bonds per domain saved as {outfiles['hbonds_prot_domains']}")


    ##########  Protein acceptors and donors per chain
    plt.figure(figsize=(10, 6))

    # Plot donor residues in red
    plt.scatter(donor_df['donor'], donor_df['frame%'], 
            color='red', label='Donor', alpha=0.6)

    # Plot acceptor residues in blue
    plt.scatter(acceptor_df['acceptor'], acceptor_df['frame%'], 
            color='blue', label='Acceptor', alpha=0.6)

    plt.xlabel('Atom Number')
    plt.ylabel('Frame %')
    plt.title('STAT1 Donor and Acceptor Residues')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig(outfiles["hbonds_prot_monomer"])
    plt.close()
    logger.info(f"Plot of protein hydrogen bonds per chain saved as {outfiles['hbonds_prot_monomer']}")

    if dna:
        ##########  dsDNA acceptors and donors 
        # Filter for DNA donors and acceptors
        dna_donor_df = hbonds_df[hbonds_df['donor_type'] == 'DNA']
        dna_acceptor_df = hbonds_df[hbonds_df['acceptor_type'] == 'DNA']

        plt.figure(figsize=(10, 6))

        # Plot DNA donor residues in red using donor atom numbers
        plt.scatter(dna_donor_df['donor'], dna_donor_df['frame%'], 
                color='red', label='Donor', alpha=0.6)

        # Plot DNA acceptor residues in blue using acceptor atom numbers
        plt.scatter(dna_acceptor_df['acceptor'], dna_acceptor_df['frame%'], 
                color='blue', label='Acceptor', alpha=0.6)

        plt.xlabel('DNA Residue Number')
        plt.ylabel('Frame %')
        plt.title('DNA Donor and Acceptor Residues')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.savefig(outfiles["hbonds_dna"])
        plt.close()
        logger.info(f"Plot of DNA hydrogen bonds saved as {outfiles['hbonds_dna']}")


        ##########  dsDNA backbone acceptors and donors
        # Filter for DNA backbone donors and acceptors
        donor_df_backbone = hbonds_df[hbonds_df['donor_section'] == 'backbone']
        acceptor_df_backbone = hbonds_df[hbonds_df['acceptor_section'] == 'backbone']

        # Create the plot
        plt.figure(figsize=(10, 6))

        # Plot donor residues in red (backbone only)
        plt.scatter(donor_df_backbone['donor'], donor_df_backbone['frame%'], 
                    color='red', label='Donor', alpha=0.6)

        # Plot acceptor residues in blue (backbone only)
        plt.scatter(acceptor_df_backbone['acceptor'], acceptor_df_backbone['frame%'], 
                    color='blue', label='Acceptor', alpha=0.6)

        # Customize the plot
        plt.xlabel('Atom Number')
        plt.ylabel('Frame %')
        plt.title('STAT1 Donor and Acceptor Residues (DNA Backbone Only)')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.savefig(outfiles["hbonds_dna_bb"])
        plt.close()
        logger.info(f"Plot of DNA backbone hydrogen bonds saved as {outfiles['hbonds_dna_bb']}")


        ##########  dsDNA bases acceptors and donors
        # Filter for DNA base donors and acceptors
        donor_df_backbone = hbonds_df[hbonds_df['donor_section'] == 'base']
        acceptor_df_backbone = hbonds_df[hbonds_df['acceptor_section'] == 'base']

        # Create the plot
        plt.figure(figsize=(10, 6))

        # Plot donor residues in red (backbone only)
        plt.scatter(donor_df_backbone['donor'], donor_df_backbone['frame%'], 
                    color='red', label='Donor', alpha=0.6)

        # Plot acceptor residues in blue (backbone only)
        plt.scatter(acceptor_df_backbone['acceptor'], acceptor_df_backbone['frame%'], 
                    color='blue', label='Acceptor', alpha=0.6)

        # Customize the plot
        plt.xlabel('Atom Number')
        plt.ylabel('Frame %')
        plt.title('STAT1 Donor and Acceptor Residues (DNA Bases Only)')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.savefig(outfiles["hbonds_dna_base"])
        plt.close()
        logger.info(f"Plot of DNA base hydrogen bonds saved as {outfiles['hbonds_dna_base']}")

if __name__ == "__main__":
    main()