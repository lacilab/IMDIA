#!/usr/bin/env python3

import sys
import pandas as pd
import argparse
from scripts.utils.helpers import setup_logger

def main(logfile: str = None, hbonds_ext: str = None, hbonds_pymol: str = None):
    """
    Generate PyMOL commands to highlight residues involved in significant hydrogen bonding.

    Arguments:
        logfile       : Path to the log file.
        hbonds_ext    : Path to the extended hydrogen bonds CSV file.
        hbonds_pymol  : Output file containing generated PyMOL commands.
    """

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    hbonds_df = pd.read_csv(hbonds_ext)

    # Filter by occupancy
    hbonds_top = hbonds_df[hbonds_df['frame%'] > 10]

    # 1. Extract unique combinations for acceptors and donors
    acceptor_combinations = hbonds_top[['acceptor_section', 'acceptor_res_num', 'acceptor_type', 'acceptor_chain']].drop_duplicates()
    donor_combinations = hbonds_top[['donor_section', 'donor_res_num', 'donor_type', 'donor_chain']].drop_duplicates()

    # Combine all unique residues (acceptors and donors might overlap)
    all_residues = pd.concat([
        acceptor_combinations.rename(columns={'acceptor_res_num': 'res_num', 'acceptor_chain': 'chain'}),
        donor_combinations.rename(columns={'donor_res_num': 'res_num', 'donor_chain': 'chain'})
    ]).drop_duplicates()


    # 2. Generate PyMOL commands
    pymol_commands = []

    # Color all atoms grey
    pymol_commands.append("color grey80, all")

    # Create selection for highlighted residues
    selection_cmd = "select bonding_residues, "
    selection_parts = []
    for _, row in all_residues.iterrows():
        res_num = int(row['res_num'])  # Ensure integer
        chain = row['chain']
        selection_parts.append(f"resi {res_num} and chain {chain}")
    selection_cmd += " or ".join(selection_parts)
    pymol_commands.append(selection_cmd)

    # Display as sticks
    pymol_commands.append("show sticks, bonding_residues")

    pymol_commands.append("color green, bonding_residues") # Color selected residues green (base color)
    pymol_commands.append('util.cnc("bonding_residues",_self=cmd)') # Color by element

    # 3. Save commands to file
    with open(hbonds_pymol, 'w') as f:
        for cmd in pymol_commands:
            f.write(cmd + '\n')
    logger.info(f"PyMOL commands saved to: {hbonds_pymol}")


if __name__ == "__main__":
    main()