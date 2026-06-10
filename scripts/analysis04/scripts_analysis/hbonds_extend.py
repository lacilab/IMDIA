#!/usr/bin/env python3

import sys
import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import MDAnalysis as mda
import argparse
from scripts.utils.helpers import setup_logger


######## HELPER FUNCTIONS ########

def get_section(res_group, res_num, domains):
    """Map atom index to domain/section."""
    if res_group not in domains:
        return 'None'
    else:
        for domain in domains[res_group]:
            if int(res_num) in range(domain['start'], domain['end']+1):
                return domain['label']
    return 'None'


def is_backbone_atom(atom_name, backbone_atoms = {'P', 'OP1', 'OP2'}):
    """Check if an atom is in the DNA backbone."""
    # Check if atom is P, OP1, OP2 or contains a prime
    return atom_name in backbone_atoms or "'" in atom_name


def get_atom_info(structure_file, atom_ids, chains):
    """Map atom index to residue number, residue name and atom name."""
    u = mda.Universe(structure_file)
    atom_info = {}
    for atom in atom_ids:
        a = u.select_atoms(f"id {atom}")
        if len(a) > 0:
            atom_info[atom] = {
                'res_name': a.resnames[0],
                'res_num': a.resids[0],
                'atom_name': a.names[0],
                'chain': get_ori_chain(a.segids[0].strip(), chains)
            }
    return atom_info

def get_ori_chain(chain_new, chains):
    """ Get original chain ID from updated ID."""
    for chain, info in chains.items():
        if info['new_id'] == chain_new:
            return chain
    return chain_new


######## MAIN ########

def main(logfile: str = None, pdbfile: str = None, universe: str = None, 
         hbonds_parsed : str = None, hbonds_ext : str = None, 
         hbonds_perframe : str = None,
         offsets : dict = None, domains : dict = {}, chaingroups : dict = {}, chains : dict = {}):

    """
    Extend hydrogen bond prediction results with data about the donor and acceptor atoms, corresponding residues and domains

    Arguments:
        logfile         : Path to the log file.
        pdbfile         : Path to input PDB file used for structural reference.
        universe        : Path to the Universe object generated with MDAnalysis, containing the processed MD system.

        hbonds_parsed   : Path to parsed hydrogen bond prediction results.
        hbonds_ext      : Path to output .csv file with extended data about donor and acceptor. Fields: 
                            donor,hydrogen,acceptor,count,frame%,
                            donor_res_name,donor_res_num,donor_atom_name,donor_chain,
                            acceptor_res_name,acceptor_res_num,acceptor_atom_name,acceptor_chain,
                            donor_type,acceptor_type,donor_res_num_adjusted,acceptor_res_num_adjusted
                            ,donor_group,acceptor_group,donor_section,acceptor_section

        hbonds_perframe : Path to output .png plot.

        offsets         : Dictionary specifying residue numbering offsets per chain.
        domains         : Dictionary of domains in the structure.
        chaingroups     : Dictionary of groups and chains in them.
        chains          : Dictionary of original ID chains and corresponding updated IDs, and offsets.
    """


    ##### Initialize logger
    logger = setup_logger("analysis", logfile)


    # Load universe
    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[0]

    # Define backbone atom criteria
    backbone_atoms = {'P', 'OP1', 'OP2'}  # Phosphate atoms without prime

    # Load parsed results
    hbonds_df = pd.read_csv(hbonds_parsed)

    # Plot hbonds per frame
    hbonds_per_frame = hbonds_df.groupby('frame').size()
    plt.figure(figsize=(10, 6))
    hbonds_per_frame.plot(kind='bar', color='skyblue')
    plt.xlabel('Frame')
    plt.ylabel('Number of Hydrogen Bonds')
    plt.title('Number of Hydrogen Bonds per Frame')
    plt.xticks(ticks=np.arange(0, len(hbonds_per_frame), 100), rotation=0)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(hbonds_perframe)
    plt.close()

    # Group by the unique bond combination and count occurrences
    hbonds_grouped = hbonds_df.groupby(['donor', 'hydrogen', 'acceptor']).size().reset_index(name='count')

    # Get total number of unique frames in the original data
    total_frames = hbonds_df['frame'].nunique()

    # Calculate percentage of frames for each bond
    hbonds_grouped['frame%'] = (hbonds_grouped['count'] / total_frames) * 100

    # Get unique atom IDs from donors and acceptors
    unique_atom_ids = set(hbonds_grouped['donor']).union(set(hbonds_grouped['acceptor']))

    # Get atom information only for relevant atoms from PDB file
    atom_info = get_atom_info(pdbfile, unique_atom_ids, chains)

    # Add donor information
    hbonds_grouped['donor_res_name'] = hbonds_grouped['donor'].map(lambda x: atom_info.get(x, {}).get('res_name', 'UNK'))
    hbonds_grouped['donor_res_num'] = hbonds_grouped['donor'].map(lambda x: atom_info.get(x, {}).get('res_num', 'UNK'))
    hbonds_grouped['donor_atom_name'] = hbonds_grouped['donor'].map(lambda x: atom_info.get(x, {}).get('atom_name', 'UNK'))
    hbonds_grouped['donor_chain'] = hbonds_grouped['donor'].map(lambda x: atom_info.get(x, {}).get('chain', 'UNK'))

    # Add acceptor information
    hbonds_grouped['acceptor_res_name'] = hbonds_grouped['acceptor'].map(lambda x: atom_info.get(x, {}).get('res_name', 'UNK'))
    hbonds_grouped['acceptor_res_num'] = hbonds_grouped['acceptor'].map(lambda x: atom_info.get(x, {}).get('res_num', 'UNK'))
    hbonds_grouped['acceptor_atom_name'] = hbonds_grouped['acceptor'].map(lambda x: atom_info.get(x, {}).get('atom_name', 'UNK'))
    hbonds_grouped['acceptor_chain'] = hbonds_grouped['acceptor'].map(lambda x: atom_info.get(x, {}).get('chain', 'UNK'))

    # Add 'donor_type' column
    hbonds_grouped['donor_type'] = np.where(
        hbonds_grouped['donor_res_name'].str.match(r'^[A-Z]{3}$'),
        'Protein',
        'DNA'
    )

    # Add 'acceptor_type' column
    hbonds_grouped['acceptor_type'] = np.where(
        hbonds_grouped['acceptor_res_name'].str.match(r'^[A-Z]{3}$'),
        'Protein',
        'DNA'
    )
    
    # Adjusted residue number to global complete protein
    offset_map = {chain: data["offset"] for chain, data in offsets.items()}
    hbonds_grouped['donor_res_num_adjusted'] = (hbonds_grouped['donor_res_num'] + hbonds_grouped['donor_chain'].map(offset_map)).astype(int)
    hbonds_grouped['acceptor_res_num_adjusted'] = (hbonds_grouped['acceptor_res_num'] + hbonds_grouped['acceptor_chain'].map(offset_map)).astype(int)

    # Chain groups
    hbonds_grouped['donor_group'] = hbonds_grouped['donor_chain'].map(lambda x: next((k for k, v in chaingroups.items() if x in v), 'None'))
    hbonds_grouped['acceptor_group'] = hbonds_grouped['acceptor_chain'].map(lambda x: next((k for k, v in chaingroups.items() if x in v), 'None'))

    # Add 'donor_section' column
    hbonds_grouped['donor_section'] = np.where(
        hbonds_grouped['donor_type'] == 'Protein',
        hbonds_grouped.apply(
        lambda row: get_section(row['donor_group'], row['donor_res_num_adjusted'], domains),
        axis=1),
        np.where(hbonds_grouped['donor_atom_name'].apply(is_backbone_atom), 'backbone', 'base')
    )

    # Add 'acceptor_section' column
    hbonds_grouped['acceptor_section'] = np.where(
        hbonds_grouped['acceptor_type'] == 'Protein',
        hbonds_grouped.apply(
        lambda row: get_section(row['acceptor_group'], row['acceptor_res_num_adjusted'], domains),
        axis=1),
        np.where(hbonds_grouped['acceptor_atom_name'].apply(is_backbone_atom), 'backbone', 'base')
    )

    # Save
    hbonds_grouped.to_csv(hbonds_ext, index=False)
    logger.info(f"Extended hydrogen bonds saved to {hbonds_ext}")


if __name__ == "__main__":
    main()