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

# Map atom index to domain/section
def get_section(res_group, res_num, domains):
    if res_group not in domains:
        return 'None'
    else:
        for domain in domains[res_group]:
            if int(res_num) in range(domain['start'], domain['end']+1):
                return domain['label']
    return 'None'


# Check if an atom is in the DNA backbone
def is_backbone_atom(atom_name, backbone_atoms = {'P', 'OP1', 'OP2'}):
    # Check if atom is P, OP1, OP2 or contains a prime
    return atom_name in backbone_atoms or "'" in atom_name


# Map atom index to residue number, residue name and atom name
def get_atom_info(structure_file, atom_ids, chains):
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
    for chain, info in chains.items():
        if info['new_id'] == chain_new:
            return chain
    return chain_new


######## MAIN ########

def main(logfile: str = None, pdbfile: str = None, universe: str = None, 
         hbonds_parsed : str = None, hbonds_ext : str = None, 
         hbonds_perframe : str = None,
         offsets : dict = None, domains : dict = {}, chaingroups : dict = {}, chains : dict = {}):

    # # Arguments
    # parser = argparse.ArgumentParser(description="Extend hydrogen bonds data with additional information.")
    # parser.add_argument("-f", "--pdb", help="PDB file.")
    # parser.add_argument("-u", "--universe", help="MDAnalysis universe file.")
    # parser.add_argument("-p", "--hbonds_parsed", help="Input. Parsed hydrogen bonds CSV file.")
    # parser.add_argument("-o", "--hbonds_ext", help="Output. Extended hydrogen bonds file.")
    # parser.add_argument("-r", "--hbonds_perframe", help="Output. Hydrogen bonds per frame file.")

    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit(0)
 
    # args = parser.parse_args()
    # pdbfile = args.pdb
    # universe = args.universe
    # hbonds_parsed = args.hbonds_parsed
    # hbonds_ext = args.hbonds_ext
    # hbonds_perframe = args.hbonds_perframe

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)


    # Load universe
    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[0]

    # Select atoms by chain
    # dimer = u.select_atoms("segid A B")  
    #monoA = u.select_atoms("segid A")   
    #monoB = u.select_atoms("segid B")   
    # dna = u.select_atoms("segid C D")   
    #dnaC = u.select_atoms("segid C") 
    #dnaD = u.select_atoms("segid D") 

    # Save indices
    # monoA_indices = set(monoA.indices)  
    # monoB_indices = set(monoB.indices)  
    # dnaC_indices = set(dnaC.indices) 
    # dnaD_indices = set(dnaD.indices) 

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