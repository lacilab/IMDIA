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
from matplotlib.ticker import MaxNLocator
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
         with_dna: bool = False, offsets: dict = None, domains: dict = {}, chaingroups: list = {}):
    """
    Quality control analysis MD results.

    Arguments:
        logfile     : Path to the log file.
        universe    : Path to the Universe object generated with MDAnalysis, containing the processed MD system.
        outdir      : Output directory for results (e.g. rep1/analysis/qc).
        outfiles    : Dictionary of the files to be generated, assigned to labels/keys.
        with_dna    : True if DNA is present in the system.
        offsets     : Dictionary specifying residue numbering offsets per chain.
        domains     : Dictionary of domains in the structure.
        chaingroups : Dictionary of groups and chains in them.
    """


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
    nucleic_sel = u.select_atoms("nucleic")
    other_sel = u.select_atoms("not (protein or nucleic)")

    ########## RMSD ##########
    rmsd_analysis = rms.RMSD(u, select='backbone', groupselections=['name CA', 'protein'])
    rmsd_analysis.run()

    rmsd_time_ns = rmsd_analysis.results.rmsd[:, 1] / 1000.0
    logger.debug(f"rmsd_analysis.results.rmsd.shape: {rmsd_analysis.results.rmsd.shape}")
    rmsd_df = pd.DataFrame(rmsd_analysis.results.rmsd[:, 2:],
                    columns=['Backbone', 'C-alphas', 'Protein'],
                    index=rmsd_time_ns)
    rmsd_df.index.name = 'Time (ns)'

    rmsd_df.to_csv(outfiles['rmsd_csv'])

    fig, ax = plt.subplots(figsize=(7, 5))
    rmsd_df.plot(title='RMSD', legend=True, ax=ax)
    ax.set_ylabel('RMSD (Å)')
    ax.set_xlabel('Time (ns)')
    ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
    ax.grid(True, alpha=0.3)
    plt.savefig(outfiles['rmsd'])
    plt.close()


    rmsd_df_bb = pd.DataFrame(
        rmsd_analysis.results.rmsd[:, 2],
        columns=['Backbone'],
        index=rmsd_time_ns
    )
    rmsd_df_bb.index.name = 'Time (ns)'

    rmsd_df_ca = pd.DataFrame(
        rmsd_analysis.results.rmsd[:, 3],
        columns=['C-alphas'],
        index=rmsd_time_ns
    )
    rmsd_df_ca.index.name = 'Time (ns)'

    rmsd_df_prot = pd.DataFrame(
        rmsd_analysis.results.rmsd[:, 4],
        columns=['Protein'],
        index=rmsd_time_ns
    )
    rmsd_df_prot.index.name = 'Time (ns)'

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    rmsd_df_bb.plot(title='RMSD Backbone', ax=axes[0], legend=False)
    rmsd_df_ca.plot(title='RMSD C-alphas', ax=axes[1], legend=False, c='orange')
    rmsd_df_prot.plot(title='RMSD Protein', ax=axes[2], legend=False, c='green')

    for ax in axes:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
        ax.set_xlabel('Time (ns)')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outfiles['rmsd_cats'])
    plt.close()

    fig, ax = plt.subplots(figsize=(7, 5))
    rmsd_df_bb.plot(title='RMSD Backbone', legend=False, ax=ax)
    ax.set_ylabel('RMSD (Å)')
    ax.set_xlabel('Time (ns)')
    ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
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
    nucleic_segids = [seg.segid for seg in nucleic_sel.segments]
    other_sesgids = [seg.segid for seg in other_sel.segments]
    logger.debug(f"Protein segids identified: {prot_segids}")
    logger.debug(f"Nucleic segids identified: {nucleic_segids}")
    logger.debug(f"Other segids identified: {other_sesgids}")

    # Calculate RMSF
    all_segids = [s for s in u.segments.segids if s.strip()] 
    logger.debug(f"All segids in the universe: {all_segids}")
    rmsf_all = {}
    for segid in all_segids:
        rmsf_all[segid] = {}
        if segid in prot_segids:
            sel = f"segid {segid} and protein and name CA"
        elif segid in nucleic_segids:
            sel = f"segid {segid} and nucleic and name P"
        else:
            logger.debug(f"Segid {segid} not recognized as protein or nucleic acid. Using backbone selection for RMSF, but this may not be ideal.")
            sel = f"segid {segid} and backbone"  # Default to all atoms in the segment, but this may not be ideal for non-protein/non-nucleic segments
        atoms = u.select_atoms(sel)
        if atoms.n_atoms == 0:
            logger.warning(f"Selection for segid {segid} is empty. Please check the chain ids.")
            continue
        else:
            logger.debug(f"Selection for segid {segid} has {atoms.n_atoms} atoms.")
            rmsf = rms.RMSF(atoms).run()
            rmsf_all[segid]['rmsf'] = rmsf.results.rmsf
            rmsf_all[segid]['resids'] = atoms.resids

    # Save protein RMSF to CSV
    rmsf_prot = pd.DataFrame()
    for chain in prot_segids:
        df_chain = pd.DataFrame({
            'Residue': rmsf_all[chain]['resids'] + offsets[chain]['offset'],   
            'Chain': chain,  
            f'RMSF': rmsf_all[chain]['rmsf']
        })
        rmsf_prot = pd.concat([rmsf_prot, df_chain], axis=1)
    
    rmsf_prot.to_csv(outfiles['rmsf_prot_csv'], index=False)

    # Plot all protein chain RMSF together
    for chain in prot_segids:
        plt.plot(rmsf_all[chain]['resids'] + offsets[chain]['offset'], rmsf_all[chain]['rmsf'], label=f'Chain {chain}')
    plt.xlabel('Residue number (adjusted)')
    plt.ylabel('RMSF (Å)')
    plt.title('RMSF of Protein Chains')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(outfiles['rmsf_prot'])
    plt.close()

    # Plot each chaingroup
    # One plot row per chaingroup, one plot per chain not in a chaingroup
    chains_with_groups = set(chain for group in chaingroups.values() for chain in group)
    chains_without_groups = [chain for chain in all_segids if chain not in chains_with_groups]
    num_groups = len(chaingroups) + len(chains_without_groups)
    logger.debug(f"Number of chaingroups: {len(chaingroups)}, Chains without groups: {len(chains_without_groups)}, Total subplots: {num_groups}")
    if num_groups > 0:
        fig, axes = plt.subplots(num_groups, 1, figsize=(10, 5*num_groups))
        if num_groups == 1:
            axes = [axes]
        for i, group in enumerate(chaingroups.keys()):
            logger.debug(f"Plotting chaingroup {group} with chains: {chaingroups[group]}")
            chain_list = chaingroups[group]
            for j, chain in enumerate(chain_list):
                if chain in rmsf_all:
                    axes[i].plot(rmsf_all[chain]['resids'] + offsets[chain]['offset'], rmsf_all[chain]['rmsf'], label=f'Chain {chain}')
                else:
                    logger.warning(f"Chain {chain} not found in rmsf_all.")
                # Color by protein domains
                if group in domains.keys():
                    group_domains = [domain for domain in domains[group]]
                    for domain in group_domains:
                        start = domain['start']
                        end = domain['end']
                        axes[i].axvspan(start, end, zorder=0, alpha=0.2, color=domain['color'], label=domain['label'])
                if j == 0:
                    axes[i].legend()
            axes[i].set_xlabel('Residue number (adjusted)')
            axes[i].set_ylabel('RMSF (Å)')
            axes[i].set_title(f'RMSF of chain group {group}')
        for i, chain in enumerate(chains_without_groups):
            j = len(chaingroups) + i
            axes[j].plot(rmsf_all[chain]['resids'] + offsets[chain]['offset'], rmsf_all[chain]['rmsf'], label=f'Chain {chain}')
            axes[j].set_xlabel('Residue number (adjusted)')
            axes[j].set_ylabel('RMSF (Å)')
            axes[j].set_title(f'RMSF of Chain {chain}')
            axes[j].legend()
        plt.tight_layout()
        plt.savefig(outfiles['rmsf_chaingroups'])
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
                'Residue': dna_rmsf[chain]['resids'] + offsets[chain]['offset'],   
                'Chain': chain,  
                f'RMSF': dna_rmsf[chain]['rmsf']
            })
            rmsf_dna = pd.concat([rmsf_dna, df_chain], axis=1)

        rmsf_dna.to_csv(outfiles['rmsf_dna_csv'], index=False)

        # Plot
        for chain in dna_rmsf.keys():
            plt.plot(dna_rmsf[chain]['resids'] + offsets[chain]['offset'], dna_rmsf[chain]['rmsf'], label=f'DNA Chain {chain}')
        plt.xlabel('Residue number')
        plt.ylabel('RMSF (Å)')
        plt.title('RMSF of DNA Chains')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(outfiles['rmsf_dna'])
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