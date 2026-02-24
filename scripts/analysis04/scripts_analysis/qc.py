#!/usr/bin/env python3

import sys
import os
import argparse
import pickle
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from MDAnalysis.analysis import rms, align
from scripts.utils.helpers import setup_logger

def radgyr(atomgroup, masses, total_mass=None):
        # coordinates change for each frame
        coordinates = atomgroup.positions
        center_of_mass = atomgroup.center_of_mass()

        # get squared distance from center
        ri_sq = (coordinates-center_of_mass)**2
        # sum the unweighted positions
        sq = np.sum(ri_sq, axis=1)
        sq_x = np.sum(ri_sq[:,[1,2]], axis=1) # sum over y and z
        sq_y = np.sum(ri_sq[:,[0,2]], axis=1) # sum over x and z
        sq_z = np.sum(ri_sq[:,[0,1]], axis=1) # sum over x and y

        # make into array
        sq_rs = np.array([sq, sq_x, sq_y, sq_z])

        # weight positions
        rog_sq = np.sum(masses*sq_rs, axis=1)/total_mass
        # square root and return
        return np.sqrt(rog_sq)

def main(logfile: str = None, universe: str = None, outdir: str = None, outfiles : dict = {}, 
         with_dna: bool = False, offsets: dict = None, domains: list = []):
    # Arguments
    # parser = argparse.ArgumentParser(description="Quality control analysis of MD simulation")
    # parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
    # parser.add_argument("-u", "--universe", help="Path to the Universe object")
    # parser.add_argument("-d", "--outdir", help="Output directory")
    # parser.add_argument("--dna", action='store_true', dest='with_dna', help="DNA is present in the system")
    # parser.add_argument("--offsets", type=dict, default=None, help="Offset for residue numbering per chain")
    # parser.add_argument("--domains", type=list, default=[], help="Protein domain definitions.")
# 
    # if len(sys.argv) < 4:
    #     parser.print_help()
    #     sys.exit(0)
# 
    # args = parser.parse_args()
    # universe = args.universe
    # outdir = args.outdir
    # with_dna = args.with_dna
    # offsets = args.offsets
    # domains = args.domains

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    os.makedirs(outdir, exist_ok=True)

    # Check if outdir exists
    if not os.path.exists(outdir):
        logger.error(f"Output directory {outdir} does not exist and could not be created.")
        sys.exit(1)

    # Check if universe file exists
    if not os.path.isfile(universe):
        logger.error(f"Universe file {universe} does not exist.")
        sys.exit(1)

    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[0]

    protein_sel = u.select_atoms("protein")

    ########## RMSD ##########
    rmsd_analysis = rms.RMSD(u, select='backbone', groupselections=['name CA', 'protein'])
    rmsd_analysis.run()
    logger.debug(f"rmsd_analysis.results.rmsd.shape: {rmsd_analysis.results.rmsd.shape}")
    rmsd_df = pd.DataFrame(rmsd_analysis.results.rmsd[:, 2:],
                    columns=['Backbone', 'C-alphas', 'Protein'],
                    index=rmsd_analysis.results.rmsd[:, 1])
    rmsd_df.index.name = 'Time (ps)'

    rmsd_df.to_csv(outfiles['rmsd_csv'])

    rmsd_df.plot(title='RMSD')
    plt.ylabel('RMSD (Å)')
    plt.savefig(outfiles['rmsd'])
    plt.close()

    rmsd_df_bb = pd.DataFrame(rmsd_analysis.results.rmsd[:, 2],
                    columns=['Backbone'],
                    index=rmsd_analysis.results.rmsd[:, 1])
    rmsd_df_bb.index.name = 'Time (ps)'

    rmsd_df_ca = pd.DataFrame(rmsd_analysis.results.rmsd[:, 3],
                        columns=['C-alphas'],
                        index=rmsd_analysis.results.rmsd[:, 1])
    rmsd_df_ca.index.name = 'Time (ps)'

    rmsd_df_prot = pd.DataFrame(rmsd_analysis.results.rmsd[:, 4],
                        columns=['Protein'],
                        index=rmsd_analysis.results.rmsd[:, 1])
    rmsd_df_prot.index.name = 'Time (ps)'

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # 1 row, 3 columns
    rmsd_df_bb.plot(title='RMSD Backbone', ax=axes[0], legend=False)
    rmsd_df_ca.plot(title='RMSD C-alphas', ax=axes[1], legend=False, c='orange')
    rmsd_df_prot.plot(title='RMSD Protein', ax=axes[2], legend=False, c='green')
    plt.tight_layout()
    plt.savefig(outfiles['rmsd_cats'])
    plt.close()

    rmsd_df_bb.plot(title='RMSD Backbone', legend=False)
    plt.ylabel('RMSD (Å)')
    plt.savefig(outfiles['rmsd_bb'])
    plt.close()

    ########## Radius of gyration ##########
    result = []
    for frame in u.trajectory:
        result.append(radgyr(atomgroup=protein_sel,
                            masses=protein_sel.masses,
                            total_mass=np.sum(protein_sel.masses)))

    result = np.asarray(result).T
    labels = ['all', 'x-axis', 'y-axis', 'z-axis']

    # Save as CSV
    df = pd.DataFrame(result, index=labels, columns=np.arange(len(u.trajectory))).T
    df.index.name = "Frame"
    df.to_csv(outfiles['radius_axis_csv'])

    for col, label in zip(result, labels):
        plt.plot(col, label=label)
    plt.legend()
    plt.ylabel('Radius of gyration (Å)')
    plt.xlabel('Frame')
    plt.title('Radius of gyration of the protein')
    plt.savefig(outfiles['radius_axis'])
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.scatter(rmsd_analysis.results.rmsd[:, 2], result[0,:],s=10, alpha=0.4, label="RMSD vs Radius of gyration")
    plt.xlabel("RMSD Backbone (Å)")
    plt.ylabel("Radius of gyration (Å)")
    plt.title("RMSD vs Radius of gyration")
    plt.legend()
    plt.savefig(outfiles['radius_rmsd'])
    plt.close()

    ########## RMSF ##########
    average = align.AverageStructure(u, u, select='protein and name CA',
                                ref_frame=0).run()
    logger.debug("Average structure calculated for RMSF analysis.")
    ref = average.results.universe

    aligner = align.AlignTraj(u, ref,
                        select='protein and name CA',
                        in_memory=True).run()
    logger.debug("Trajectory aligned to average structure for RMSF analysis.")

    # Get segids of the protein, to avoid leaving out non-standard residues like PTR
    prot_segids = [seg.segid for seg in protein_sel.segments]
    logger.debug(f"Protein segids identified: {prot_segids}")

    prot_rmsf = {}
    for segid in prot_segids:
        prot_rmsf[segid] = {}
        ca = u.select_atoms(f"segid {segid} and name CA")
        prot_rmsf[segid]['resids'] = ca.resids
        if ca.n_atoms == 0:
            logger.warning(f"Selection for segid {segid} is empty. Please check the chain ids.")
            continue
        else:
            logger.debug(f"Selection for segid {segid} has {ca.n_atoms} CA atoms.")
            rmsf = rms.RMSF(ca).run()
            prot_rmsf[segid]['rmsf'] = rmsf.results.rmsf

    # Save Protein RMSF to CSV
    rmsf_prot = pd.DataFrame()
    for chain in prot_rmsf.keys():  
        df_chain = pd.DataFrame({
            'Residue': prot_rmsf[chain]['resids'] + offsets[chain],   
            'Chain': chain,  
            f'RMSF': prot_rmsf[chain]['rmsf']
        })
        rmsf_prot = pd.concat([rmsf_prot, df_chain], axis=1)
    
    rmsf_prot.to_csv(outfiles['rmsf_prot_csv'], index=False)

    # Plot
    for chain in prot_rmsf.keys():
        plt.plot(prot_rmsf[chain]['resids'] + offsets[chain], prot_rmsf[chain]['rmsf'], label=f'Chain {chain}')
    # Color by protein domains
    if len(domains) > 0:
        for domain in domains:
            start = domain['start']
            end = domain['end']
            plt.axvspan(start, end, zorder=0, alpha=0.2, color=domain['color'], label=domain['label'])
    plt.xlabel('Residue number (adjusted)')
    plt.ylabel('RMSF (Å)')
    plt.title('RMSF of Protein Chains')
    plt.legend()
    plt.savefig(outfiles['rmsf_prot_domains'])
    plt.close()

    ##### If the complex contains dsDNA #####
    if with_dna:
        dna_sel = u.select_atoms("nucleic") 
        dna_segids = [seg.segid for seg in dna_sel.segments]

        dna_rmsf = {}
        for segid in dna_segids:
            dna_rmsf[segid] = {}
            p = u.select_atoms(f"segid {segid} and name P")
            dna_rmsf[segid]['resids'] = p.resids
            if p.n_atoms == 0:
                logger.warning(f"Selection for segid {segid} is empty. Please check the chain ids.")
                continue
            else:
                logger.debug(f"Selection for segid {segid} has {p.n_atoms} P atoms.")
                rmsf = rms.RMSF(p).run()
                dna_rmsf[segid]['rmsf'] = rmsf.results.rmsf

        # Save DNA RMSF to CSV
        rmsf_dna = pd.DataFrame()
        for chain in dna_rmsf.keys():  
            df_chain = pd.DataFrame({
                'Residue': dna_rmsf[chain]['resids'] + offsets[chain],   
                'Chain': chain,  
                f'RMSF': dna_rmsf[chain]['rmsf']
            })
            rmsf_dna = pd.concat([rmsf_dna, df_chain], axis=1)

        rmsf_dna.to_csv(outfiles['rmsf_dna_csv'], index=False)

        # Plot
        for chain in dna_rmsf.keys():
            plt.plot(dna_rmsf[chain]['resids'] + offsets[chain], dna_rmsf[chain]['rmsf'], label=f'DNA Chain {chain}')
        plt.xlabel('Residue number')
        plt.ylabel('RMSF (Å)')
        plt.title('RMSF of DNA Chains')
        plt.legend()
        plt.savefig(outfiles['rmsf_dna_res'])
        plt.close()

    logger.debug("\nChecking generated files.")

    expected_files = outfiles.values()

    missing = [p for p in expected_files if not os.path.isfile(p)]
    for file in expected_files:
        if file not in missing:
            logger.debug(f"generated {file}")
        else:
            logger.warning(f"{file} not generated")

    if not missing:
        logger.info("QC completed successfully. All expected files were generated.")
    else:
        logger.warning("QC completed with warnings. Some expected files are missing.")

if __name__ == "__main__":
    main()