#!/usr/bin/env python3

import pandas as pd
from scripts.utils.helpers import setup_logger

def main(logfile: str = None, hbonds_gmx: str = None, hbonds_gmx_parsed: str = None, group1: str = None, group2: str = None):
    # # Arguments
    # parser = argparse.ArgumentParser(description="Parse GROMACS hydrogen bonds output into a CSV file.")
    # parser.add_argument("-i", "--ingmx", help="Input file, from the GROMACS hbond analysis.")
    # parser.add_argument("-o", "--outcsv", help="Output path to save the parsed CSV file.")
# 
    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit(0)
# 
    # args = parser.parse_args()
    # hbonds_gmx = args.ingmx
    # hbonds_gmx_parsed = args.outcsv

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