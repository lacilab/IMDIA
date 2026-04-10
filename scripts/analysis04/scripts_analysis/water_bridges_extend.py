#!/usr/bin/env python3

from pathlib import Path
import pickle
import pandas as pd
import numpy as np
from scripts.utils.helpers import setup_logger


def get_original_chain_id(universe, atom_number, chain_dict):
    """Full pipeline: atom → chainID → original chain key."""
    # 1. MDAnalysis: atom number → chainID
    atom = universe.atoms[atom_number - 1]
    new_id = atom.chainID  
    # 2. Dict lookup: new_id → original chain key
    if len(new_id.strip()) > 0:
        old_id = next((key for key, val in chain_dict.items() if val['new_id'] == new_id), None)
    else:
        old_id = 'UNK'

    return old_id


def main(logfile: str, universe: str, wb_csv: Path, offsets: dict, outfile: Path):

    # Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug("Parsing and filtering water bridge analysis results.")

    df = pd.read_csv(wb_csv, sep=',')

    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[0]

    # Assign chain IDs
    df["sele1_chain"] = df['sele1_index'].apply(lambda atom_num: get_original_chain_id(u, atom_num, offsets))
    df["sele2_chain"] = df['sele2_index'].apply(lambda atom_num: get_original_chain_id(u, atom_num, offsets))

    # Adjusted residue numbers based on offsets
    df["sele1_resid_adj"] = (
        df['sele1_resid'] + 
        df['sele1_chain'].map(lambda chain: offsets[chain]['offset'] if chain != 'UNK' and pd.notna(chain) else 0)
    )
    df["sele2_resid_adj"] = (
        df['sele2_resid'] + 
        df['sele2_chain'].map(lambda chain: offsets[chain]['offset'] if chain != 'UNK' and pd.notna(chain) else 0)
    )

    # Filter out entries where none of the residues is WAT
    df = df[~((df['sele1_resnm'] != 'WAT') & (df['sele2_resnm'] != 'WAT'))]

    # Aggregate by atom pair
    agg_df = df.groupby(['sele1_index', 'sele2_index', 'sele1_chain', 'sele2_chain', 'sele1_resid', 'sele2_resid', 'sele1_resid_adj', 
                          'sele2_resid_adj', 'sele1_resnm', 'sele2_resnm']).size().reset_index(name='count')

    # Save extended df to csv
    agg_df.to_csv(outfile, index=False)



if __name__ == "__main__":
    main()