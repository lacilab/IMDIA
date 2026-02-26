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
def get_section(res_num, prot_domains):
    if pd.isna(res_num): 
        return 'Unknown' 
    for domain, res_range in prot_domains.items():
        if res_num in res_range:
            return domain
    return 'Unknown_protein'


# Check if an atom is in the DNA backbone
def is_backbone_atom(atom_name, backbone_atoms = {'P', 'OP1', 'OP2'}):
    # Check if atom is P, OP1, OP2 or contains a prime
    return atom_name in backbone_atoms or "'" in atom_name


# Map atom index to residue number, residue name and atom name
def get_atom_info(structure_file, atom_ids, chains=['A', 'B', 'C', 'D']):
    atom_info = {}
    with open(structure_file, 'r') as f:
        for line in f:
            # Check for ATOM or HETATM records
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    atom_id = int(line[6:11].strip())      # Serial number (columns 7-11)
                    chain= line[21].strip()            # Chain ID (column 22)
                    
                    if atom_id in atom_ids and chain in chains:
                        atom_name = line[12:16].strip()    # Atom name (columns 13-16)
                        res_name = line[17:20].strip()     # Residue name (columns 18-20)
                        res_num = int(line[22:26].strip()) # Residue sequence number (columns 23-26)                        
                        # Store in dictionary
                        atom_info[atom_id] = {
                            'res_name': res_name,
                            'res_num': res_num,
                            'atom_name': atom_name,
                            'chain': chain
                        }
                except (ValueError, IndexError):
                    continue
    
    return atom_info

######## MAIN ########

def main(logfile: str = None, pdbfile: str = None, universe: str = None, 
         hbonds_parsed : str = None, hbonds_ext : str = None, 
         hbonds_perframe : str = None,
         offsets : dict = None, domains : list = []):

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

    # Protein domains
    prot_domains = {}
    if len(domains) > 0:
        for domain in domains:
            prot_domains[domain['label']] = range(domain['start'], domain['end'])


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
    atom_info = get_atom_info(pdbfile, unique_atom_ids)

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
        'protein',
        'DNA'
    )

    # Add 'acceptor_type' column
    hbonds_grouped['acceptor_type'] = np.where(
        hbonds_grouped['acceptor_res_name'].str.match(r'^[A-Z]{3}$'),
        'protein',
        'DNA'
    )
    
    # Adjusted residue number to global complete protein
    offset_map = {chain: data["offset"] for chain, data in offsets.items()}
    hbonds_grouped['donor_res_num_adjusted'] = (hbonds_grouped['donor_res_num'] + hbonds_grouped['donor_chain'].map(offset_map)).astype(int)
    hbonds_grouped['acceptor_res_num_adjusted'] = (hbonds_grouped['acceptor_res_num'] + hbonds_grouped['acceptor_chain'].map(offset_map)).astype(int)

    # Add 'donor_section' column
    hbonds_grouped['donor_section'] = np.where(
        hbonds_grouped['donor_type'] == 'protein',
        hbonds_grouped.apply(
        lambda row: get_section(row['donor_res_num_adjusted'], prot_domains),
        axis=1),
        np.where(hbonds_grouped['donor_atom_name'].apply(is_backbone_atom), 'backbone', 'base')
    )

    # Add 'acceptor_section' column
    hbonds_grouped['acceptor_section'] = np.where(
        hbonds_grouped['acceptor_type'] == 'protein',
        hbonds_grouped.apply(
        lambda row: get_section(row['acceptor_res_num_adjusted'], prot_domains),
        axis=1),
        np.where(hbonds_grouped['acceptor_atom_name'].apply(is_backbone_atom), 'backbone', 'base')
    )

    # Save
    hbonds_grouped.to_csv(hbonds_ext, index=False)
    logger.info(f"Extended hydrogen bonds saved to {hbonds_ext}")


if __name__ == "__main__":
    main()