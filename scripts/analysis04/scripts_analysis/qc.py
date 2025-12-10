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

def main(logfile: str = None, universe: str = None, outdir: str = None, with_dna: bool = False, offset: int = 0, domains: list = []):
    # Arguments
    # parser = argparse.ArgumentParser(description="Quality control analysis of MD simulation")
    # parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
    # parser.add_argument("-u", "--universe", help="Path to the Universe object")
    # parser.add_argument("-d", "--outdir", help="Output directory")
    # parser.add_argument("--dna", action='store_true', dest='with_dna', help="DNA is present in the system")
    # parser.add_argument("--offset", type=int, default=0, help="Offset for residue numbering")
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
    # offset = args.offset
    # domains = args.domains

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    # Output file paths
    rmsd = f"{outdir}/rmsd.png"
    rmsd_cats = f"{outdir}/rmsd_cats.png"
    rmsd_bb = f"{outdir}/rmsd_bb.png"
    radius_axis = f"{outdir}/radius_axis.png"
    radius_rmsd = f"{outdir}/radius_rmsd.png"
    rmsf_prot_domains = f"{outdir}/rmsf_prot_domains.png"
    rmsf_dna_res = f"{outdir}/rmsf_dna_res.png"
    rmsd_csv = f"{outdir}/rmsd.csv"
    radius_axis_csv = f"{outdir}/radius_axis.csv"
    rmsf_prot_csv = f"{outdir}/rmsf_prot.csv"
    rmsf_dna_csv = f"{outdir}/rmsf_dna.csv"

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

    # Select atoms by chain
    dimer = u.select_atoms("segid A B")  
    monoA = u.select_atoms("segid A")   
    monoB = u.select_atoms("segid B")

    # Check if protein selections are empty
    prot_sele = {"dimer": dimer, "monoA": monoA, "monoB": monoB}
    for name, sele in prot_sele.items():
        if sele.n_atoms == 0:
            logger.error(f"Selection {name} is empty. Please check the chain ids.")
            sys.exit(1)

    ########## RMSD ##########
    rmsd_analysis = rms.RMSD(u, select='backbone', groupselections=['name CA', 'protein'])
    rmsd_analysis.run()
    logger.debug(f"rmsd_analysis.results.rmsd.shape: {rmsd_analysis.results.rmsd.shape}")
    rmsd_df = pd.DataFrame(rmsd_analysis.results.rmsd[:, 2:],
                    columns=['Backbone', 'C-alphas', 'Protein'],
                    index=rmsd_analysis.results.rmsd[:, 1])
    rmsd_df.index.name = 'Time (ps)'

    rmsd_df.to_csv(rmsd_csv)

    rmsd_df.plot(title='RMSD')
    plt.ylabel('RMSD (Å)')
    plt.savefig(rmsd)
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
    plt.savefig(rmsd_cats)
    plt.close()

    rmsd_df_bb.plot(title='RMSD Backbone', legend=False)
    plt.ylabel('RMSD (Å)')
    plt.savefig(rmsd_bb)
    plt.close()

    ########## Radius of gyration ##########
    

    result = []
    for frame in u.trajectory:
        result.append(radgyr(atomgroup=dimer,
                            masses=dimer.masses,
                            total_mass=np.sum(dimer.masses)))

    result = np.asarray(result).T
    labels = ['all', 'x-axis', 'y-axis', 'z-axis']

    # Save as CSV
    df = pd.DataFrame(result, index=labels, columns=np.arange(len(u.trajectory))).T
    df.index.name = "Frame"
    df.to_csv(radius_axis_csv)

    for col, label in zip(result, labels):
        plt.plot(col, label=label)
    plt.legend()
    plt.ylabel('Radius of gyration (Å)')
    plt.xlabel('Frame')
    plt.savefig(radius_axis)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.scatter(rmsd_analysis.results.rmsd[:, 2], result[0,:],s=10, alpha=0.4, label="RMSD vs Radius of gyration")
    plt.xlabel("RMSD Backbone (Å)")
    plt.ylabel("Radius of gyration (Å)")
    plt.title("RMSD vs Radius of gyration")
    plt.legend()
    plt.savefig(radius_rmsd)
    plt.close()

    ########## RMSF ##########
    average = align.AverageStructure(u, u, select='protein and name CA',
                                ref_frame=0).run()
    ref = average.results.universe

    aligner = align.AlignTraj(u, ref,
                        select='protein and name CA',
                        in_memory=True).run()

    monoA_ca = monoA.select_atoms('name CA')
    RA = rms.RMSF(monoA_ca).run()
    monoB_ca = monoB.select_atoms('name CA')
    RB = rms.RMSF(monoB_ca).run()

    # Save Protein RMSF to CSV
    rmsf_prot = pd.DataFrame({
        'Residue': monoA_ca.resids,
        'RMSF_A': RA.results.rmsf,
        'RMSF_B': RB.results.rmsf
    })
    rmsf_prot.to_csv(rmsf_prot_csv, index=False)

    # Plot
    plt.plot(monoA_ca.resids + offset, RA.results.rmsf, label='RMSF Chain A')
    plt.plot(monoB_ca.resids + offset, RB.results.rmsf, label='RMSF Chain B')
    # Color by protein domains
    if len(domains) > 0:
        for domain in domains:
            start = domain['start']
            end = domain['end']
            plt.axvspan(start, end, zorder=0, alpha=0.2, color=domain['color'], label=domain['label'])
    plt.xlabel('Residue number (adjusted)')
    plt.ylabel('RMSF (Å)')
    plt.legend()
    plt.savefig(rmsf_prot_domains)
    plt.close()

    ##### If the complex contains dsDNA #####
    if with_dna:
        dna = u.select_atoms("segid C D")   
        dnaC = u.select_atoms("segid C") 
        dnaD = u.select_atoms("segid D") 

        # Select P atoms to focus on backbone fluctuations
        dnaC_p = dnaC.select_atoms('name P')
        RC = rms.RMSF(dnaC_p).run()
        dnaD_p = dnaD.select_atoms('name P')
        RD = rms.RMSF(dnaD_p).run()

        # Check if DNA selections are empty
        dna_sele = {"dnaC_p": dnaC_p, "dnaD_p": dnaD_p}
        for name, sele in dna_sele.items():
            if sele.n_atoms == 0:
                logger.error(f"Selection {name} is empty. Please check the chain ids.")
                sys.exit(1)

        # Save DNA RMSF to CSV
        rmsf_dna = pd.DataFrame({
            'Residue': dnaC_p.resids,
            'RMSF_C': RC.results.rmsf,
            'RMSF_D': RD.results.rmsf
        })
        rmsf_dna.to_csv(rmsf_dna_csv, index=False)

        # Plot
        plt.plot(dnaC_p.resids, RC.results.rmsf, label='RMSF DNA C')
        plt.plot(dnaD_p.resids, RD.results.rmsf, label='RMSF DNA D')
        plt.xlabel('Residue number')
        plt.ylabel('RMSF (Å)')
        plt.legend()
        plt.savefig(rmsf_dna_res)
        plt.close()

    logger.debug("\nChecking generated files.")

    expected_files_prot = [
        rmsd,
        rmsd_cats,
        rmsd_bb,
        radius_axis,
        radius_rmsd,
        rmsf_prot_domains,
        rmsd_csv,
        radius_axis_csv,
        rmsf_prot_csv
    ]
    expected_files_dna = [
        rmsf_dna_res,
        rmsf_dna_csv
    ]

    expected_files = (expected_files_prot + expected_files_dna, expected_files_prot)[with_dna]

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