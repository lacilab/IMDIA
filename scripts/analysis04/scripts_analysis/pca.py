#!/usr/bin/env python3

import sys
import os
import pickle
import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import pca, align
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scripts.utils.helpers import setup_logger


def main(logfile : str = None, universe: str = None, outdir: str = None):
    # Arguments
    # parser = argparse.ArgumentParser(description="Analyze water bridges between protein and nucleic acids.")
    # parser.add_argument("-u", "--universe", help="Path to the Universe object")
    # parser.add_argument("-d", "--outdir", help="Directory for output data")
# 
    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit(0)
# 
    # args = parser.parse_args()
    # outdir = args.outdir
    # universe = args.universe

    ##### Initialize logger
    logger = setup_logger("analysis", logfile)

    # Output file paths
    pc_cumulative = f"{outdir}/pc_cumulative.png"
    pc_pairgrid = f"{outdir}/pc_pairgrid.png"
    pc_cosine = f"{outdir}/pc_cosine.png"

    os.makedirs(outdir, exist_ok=True)

    with open(universe, 'rb') as f:
        u = pickle.load(f)
    u.trajectory[0]

    aligner = align.AlignTraj(u, u, select='protein and name CA',
                            in_memory=True).run()

    pc = pca.PCA(u, select='protein and name CA',
                align=True, mean=None,
                n_components=None).run()

    backbone = u.select_atoms('protein and name CA')
    n_bb = len(backbone)
    logger.debug(f'There are {n_bb} backbone atoms in the analysis')
    logger.debug(f'PC components shape: {pc.results.p_components.shape}')

    # Plot cumulative variance captured by PCs
    plt.plot(pc.results.cumulated_variance[:15])
    plt.xlabel('Principal component')
    plt.ylabel('Cumulative variance')
    plt.savefig(pc_cumulative)
    plt.close()
    logger.info(f"Cumulative variance plot saved to: {pc_cumulative}")

    # Plot pairgrid
    transformed = pc.transform(backbone, n_components=3)
    df_pca = pd.DataFrame(transformed,
                    columns=['PC{}'.format(i+1) for i in range(3)])
    df_pca['Time (ps)'] = df_pca.index * u.trajectory.dt
    g = sns.PairGrid(df_pca, hue='Time (ps)',
                    palette=sns.color_palette('Oranges_d',
                                            n_colors=len(df_pca)))
    g.map(plt.scatter, marker='.')
    plt.savefig(pc_pairgrid)
    plt.close()
    logger.info(f"Pairgrid plot saved to: {pc_pairgrid}")

    # Cosine content
    # 0 to 1, <0.7 is OK
    for i in range(3):
        cc = pca.cosine_content(transformed, i)
        logger.debug(f"Cosine content for PC {i+1} = {cc:.3f}")

    melted = pd.melt(df_pca, id_vars=["Time (ps)"],
                    var_name="PC",
                    value_name="Value")
    g = sns.FacetGrid(melted, col="PC")
    g.map(sns.lineplot,
        "Time (ps)", # x-axis
        "Value", # y-axis
        errorbar=None) # no confidence interval
    plt.savefig(pc_cosine)
    plt.close()
    logger.info(f"Cosine content plot saved to: {pc_cosine}")

if __name__ == "__main__":
    main()
