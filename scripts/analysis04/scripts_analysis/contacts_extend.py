#!/usr/bin/env python3

from pathlib import Path
import pickle
import pandas as pd
import numpy as np
from scripts.utils.helpers import setup_logger


def main(logfile: str, universe: str, contacts_csv: Path, offsets: dict, outfile: Path):
    # Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug("Parsing contact analysis results.")

    # Load universe
    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[-1] 
    nframes = u.trajectory.frame

    df = pd.read_csv(contacts_csv, sep=',')

 
    # Adjusted residue numbers based on offsets
    df["sele1_resid_adj"] = (
        df['sele1_resid'] + 
        df['sele1_chain'].map(lambda chain: offsets[chain]['offset'] if chain in offsets.keys() else 0)
    )
    df["sele2_resid_adj"] = (
        df['sele2_resid'] + 
        df['sele2_chain'].map(lambda chain: offsets[chain]['offset'] if chain in offsets.keys() else 0)
    )

    # Aggregate by residue pair, with avg distance and count of interactions / number of frames 
    agg_df = df.groupby(['sele1_resnm', 'sele1_resid_adj', 'sele1_chain', 
                              'sele2_resnm', 'sele2_resid_adj', 'sele2_chain']).agg(
        count=('distance', 'size'), 
        distance_avg=('distance', 'mean') 
    ).reset_index()

    agg_df['ratio'] = agg_df['count'] / nframes

    # Save extended df to csv
    agg_df.to_csv(outfile, index=False)


if __name__ == "__main__":
    main()