#!/usr/bin/env python3

import sys
import MDAnalysis as mda
import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')
import os
import pickle
# import argparse
from scripts.utils.helpers import setup_logger

def main(logfile: str = None, pdbfile : str = None, xtcfile : str = None, outdir : str = None, universe : str = None):
    """
    Create an MDAnalysis Universe from PDB and XTC files for subsequent trajectory analysis.

    Arguments:
        logfile   : Path to the log file.
        pdbfile   : Path to the input PDB structure file.
        xtcfile   : Path to the MD trajectory file (.xtc).
        outdir    : Output directory for generated files (e.g. rep1/analysis/raw).
        universe  : Path to save MDAnalysis Universe object.
    """

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    os.makedirs(outdir, exist_ok=True)

    u = mda.Universe(pdbfile, xtcfile)
    with open(universe, 'wb') as f:
        pickle.dump(u, f)

    logger.debug(f'Number of frames in the trajectory: {len(u.trajectory)}')
    logger.debug(f'Number of residues: {len(u.residues)}')
    logger.debug(f'Number of atoms: {len(u.atoms)}')
    logger.debug(f'Segments: {u.segments}')

if __name__ == "__main__":
    main()