#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
from scripts.utils.helpers import setup_logger


def main(logfile, hbonds_files, outfile):
    """ 
    Merge hydrogen bond prediction results from different tools / .csv files. 
    
    Arguments:
        logfile       : Path to the log file.
        hbonds_files  : List of .csv files to merge.
        outfile       : Output .csv file.
    """

    
    # Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug("Merging hbond extended files.")

    # Read and concatenate all hbond files
    all_hbonds = []
    for file in hbonds_files:
        logger.debug(f"Reading hbond file: {file}")
        # Check file exists and is .csv
        if Path(file).is_file() and Path(file).suffix == ".csv":
            df = pd.read_csv(file)
            all_hbonds.append(df)
    merged_hbonds = pd.concat(all_hbonds, axis=0)
    logger.debug(f"Merged {len(hbonds_files)} hbond files into a single DataFrame with {len(merged_hbonds)} entries")


    # Save and check
    try:
        merged_hbonds.to_csv(outfile, index=False)
        logger.info(f"Merged hydrogen bond analysis results saved to {outfile}")      
    except OSError as err:
        logger.warning(f"Merged hydrogen bond analysis results could not be saved. Check log {logfile}")
        logger.warning(f"Check log {logfile}")
    return 

if __name__ == "__main__":
    main()