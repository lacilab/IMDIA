#!/usr/bin/env python3

import sys
import os
import re
import pickle
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import numpy as np
import pandas as pd
from scripts.utils.helpers import setup_logger
from scripts.utils.parallel import an_perblock

def hbond_factory(u, start, stop, sele_index, dha_angle, da_angle):
    hbonds = HBA(universe=u, between=[sele_index], 
                    d_h_a_angle_cutoff=dha_angle, d_a_cutoff=da_angle)
    hbonds.run(start=start, stop=stop, verbose=True)
    return hbonds.results.hbonds


def main(logfile: str = None, universe: str = None, 
         tprfile: str = None, xtcfile: str = None, 
         outdir : str = None, hbonds_mda_parsed : str = None,
         cores : int = 1,
         start: int = 0, stop: int = None,
         group1: str = None, group2: str = None, dna : bool = False, offsets : dict = None):

    """
    Predict hydrogen bonds between two groups using MDAnalysis.

    Arguments:
        logfile           : Path to the log file.
        universe          : Path to the Universe object generated with MDAnalysis, containing the processed MD system.
        tprfile           : GROMACS run input file (.tpr) containing topology and run parameters.
        xtcfile           : GROMACS trajectory file (.xtc) containing atomic coordinates over time.
        outdir            : Output directory for results (e.g. rep1/analysis/hbonds).
        hbonds_mda_parsed : Output .csv with the parsed results, with fields frame,donor,hydrogen,acceptor,distance,angle.

        cores             : Number of CPU cores to use for parallel processing.

        start             : Start frame for analysis window.
        stop              : End frame for analysis window.
       
        group1            : First atom group for hydrogen bond analysis.
        group2            : Second atom group for hydrogen bond analysis.

        dna               : True if DNA is present in the system.
        offsets           : Dictionary specifying residue numbering offsets per chain.
    """

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug("Starting hydrogen bond analysis with MDAnalysis.")

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
        logger.info(f"Created output directory: {outdir}")
    else:
        logger.warning(f"Output directory already exists: {outdir}")

    # Check that group1 and group2 are provided
    if group1 is None or group2 is None:
        logger.error("Both group1 and group2 must be specified for interaction analysis.")
        sys.exit(1)

    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[0]

    u_tpr = mda.Universe(tprfile,xtcfile)

    # Define selections
    selections = []
    for group in [group1, group2]:
        if re.match(r"^[A-Za-z]$", group):
            new_id = next((k for k, v in offsets.items() if v['id_old'] == group), 'UNK')
            logger.debug(f"MDA - hbond analysis - Group {group} identified as chain, with new id {new_id}.")
            selections.append(u.select_atoms(f"segid {new_id}"))
        elif group == "DNA" and dna:
            logger.debug(f"MDA - hbond analysis - Group {group} identified as DNA.")
            selections.append(u.select_atoms("nucleic")) 
        elif group == "Protein":
            logger.debug(f"MDA - hbond analysis - Group {group} identified as Protein.")
            protein_sel = u.select_atoms("protein")
            # Get segids of the protein, to avoid leaving out non-standard residues like PTR
            prot_segids = [seg.segid for seg in protein_sel.segments]
            logger.debug(f"Protein segids identified: {prot_segids}")
            selections.append(u.select_atoms("segid " + " ".join(prot_segids)))

    # Convert AtomGroup indices to selection strings
    sele_index = []
    for sel in selections:
        sele_index.append("index " + " ".join(map(str, sel.indices)))

    # Define angle cutoffs
    dha_angle, da_angle = 120, 3.5


    if cores > 1:
        logger.debug(f"MDA - hbond analysis - Running parallel Dask with {cores} cores")
        hbonds_array = an_perblock(logger, tprfile, xtcfile, cores, hbond_factory, start, stop, sele_index, dha_angle, da_angle)
    else:
        logger.debug(f"MDA - hbond analysis - Running single-core analysis.")
        hbonds = HBA(
            universe=u_tpr,
            between=[sele_index],
            d_h_a_angle_cutoff=dha_angle,
            d_a_cutoff=da_angle
        )
        # Run analysis
        hbonds.run(start=start, stop=stop, verbose=True)
        # Extract results
        hbonds_array = hbonds.results.hbonds

    # Name columns and save raw data
    df = pd.DataFrame(hbonds_array, columns=["frame", "donor", "hydrogen", "acceptor", "distance", "angle"])
    # Adjust to 1-index
    df['donor'] = (df['donor'] + 1).astype('int')
    df['hydrogen'] = (df['hydrogen'] + 1).astype('int')
    df['acceptor'] = (df['acceptor'] + 1).astype('int')

    # Save
    df.to_csv(hbonds_mda_parsed, index=False)
    logger.info(f"Parsed hydrogen bonds saved to {hbonds_mda_parsed}")


if __name__ == "__main__":
    main()