#!/usr/bin/env python3

import pandas as pd
from scripts.utils.helpers import setup_logger

def main(logfile: str = None, hbonds_gmx: str = None, hbonds_gmx_parsed: str = None, group1: str = None, group2: str = None):
    """
    Predict hydrogen bonds between two groups using GROMACS analysis tools.

    Arguments:
        logfile           : Path to the log file.
        hbonds_gmx        : Path to GROMACS .ndx file generated during the hbond analysis (e.g. analysis/hbonds/hbonds_gmx.ndx)
        hbonds_gmx_parsed : Output .csv file with donor,hydrogen,acceptor,frame fields.
        group1            : First group used in hydrogen bond analysis.
        group2            : Second group used in hydrogen bond analysis.
    """

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    # Check that group1 and group2 are provided
    if group1 is None or group2 is None:
        logger.error("Both group1 and group2 must be specified to parse the gmx hbonds output.")
        sys.exit(1)

    # Lists to store hydrogen bond data
    hbonds_data = []
    current_frame = None

    # Parse the .ndx file
    with open(hbonds_gmx, 'r') as f:
        lines = f.readlines()
        start_parsing = False
        
    for line in lines:
        line = line.strip()
        
        # Start parsing at [ hbonds_Chain_A-Chain_B_frame_0 ]
        if f'[ hbonds_{group1}-{group2}_frame_0 ]' in line:
            start_parsing = True
            current_frame = 0
            continue
        
        # Stop parsing if hit an unrelated section
        if start_parsing and line.startswith('[') and f'hbonds_{group1}-{group2}_frame' not in line:
            break
        
        if start_parsing:
            # Check if line is a new frame header
            if line.startswith(f'[ hbonds_{group1}-{group2}_frame_'):
                frame_num = int(line.split('_')[-1].replace(' ]', ''))
                current_frame = frame_num
            # Process lines with atom indices
            elif line and not line.startswith('['):
                try:
                    donor, hydrogen, acceptor = map(int, line.split())
                    hbonds_data.append({
                        'donor': donor,
                        'hydrogen': hydrogen,
                        'acceptor': acceptor,
                        'frame': current_frame
                    })
                except ValueError:
                    continue  # Skip malformed lines

    df = pd.DataFrame(hbonds_data)
    df.to_csv(hbonds_gmx_parsed, index=False)
    logger.info(f"Parsed hydrogen bonds saved to {hbonds_gmx_parsed}")

if __name__ == "__main__":
    main()