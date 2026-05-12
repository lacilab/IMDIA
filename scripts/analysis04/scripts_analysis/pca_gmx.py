from pathlib import Path
import sys
import subprocess
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scripts.utils.helpers import setup_logger

def main(logfile: Path, outdir: Path, tpr: Path, xtc: Path):
   """
   Run principal component analysis (PCA) on MD results using GROMACS.

   Arguments:
      logfile  : Path to the log file.
      outdir   : Output directory for results (e.g. rep1/analysis/pca).
      tpr  : GROMACS run input file (.tpr) containing topology and run parameters.
      xtc  : GROMACS trajectory file (.xtc) containing atomic coordinates over time.
   """

   if not logfile.is_file():
      print(f"Log file {logfile} does not exist.")
      sys.exit(1)

   logger = setup_logger("analysis", logfile)

   outdir.mkdir(parents=True, exist_ok=True)

   if not tpr.is_file():
      logger.error(f"TPR file {tpr} does not exist.")
      sys.exit(1)
   if not xtc.is_file():
      logger.error(f"XTC file {xtc} does not exist.")
      sys.exit(1)

   # Build and diagonalize the covariance matrix
   cmd = ["gmx", "covar", "-f", str(xtc), "-s", str(tpr), "-o", str(outdir / "eigenval.xvg"), 
          "-v", str(outdir / "eigenvec.trr"), "-av", str(outdir / "average.pdb"), "-l", str(outdir / "covar.log"), "-quiet"]
   subprocess.run(cmd, check=True, text=True, input="C-alpha\nC-alpha\n")

   # Analyze cosine content
   cmd = ["gmx", "anaeig", "-f", str(xtc), 
          "-v", str(outdir / "eigenvec.trr"), "-eig", str(outdir / "eigenval.xvg"), "-s", str(tpr), 
          "-proj", str(outdir / "proj_10.xvg"), "-rmsf", str(outdir / "eigrmsf.xvg"), "-b", "0", "-first", "1", "-last", "10", "-quiet"]
   subprocess.run(cmd, check=True, text=True, input="C-alpha\nC-alpha\n")
   cmd = ["gmx", "analyze", "-f", str(outdir / "proj_10.xvg"), "-n", "10", "-cc", str(outdir / "pc_cosine_corr.xvg"), "-quiet"]
   subprocess.run(cmd, check=True, text=True)
   
   cc = np.loadtxt(outdir / "pc_cosine_corr.xvg", comments=['@','#'])
   plt.plot(cc[:,0], cc[:,1], 'o-')
   plt.xlabel('Principal component')
   plt.gca().set_xticks(np.arange(1, 11, 1)) 
   plt.ylabel('Cosine content')
   plt.grid(True, alpha=0.3)
   plt.title('Cosine content of PCs')
   plt.savefig(outdir / 'pc_cosine.png')
   plt.close()

   # Plot cummulative variance captured by PCs
   ev = np.loadtxt(outdir /'eigenval.xvg', comments=['@','#'], usecols=[1])
   total_var = ev.sum()
   cum_var = np.cumsum(ev) / total_var * 100
   pc_lim = min(20, len(ev)) # Plot first 20 PCs
   plt.plot(range(1, len(ev) + 1), cum_var, 'o-')
   plt.xlabel('Principal component')
   plt.xlim(0, pc_lim + 1)
   plt.gca().set_xticks(np.arange(1, pc_lim + 1, 1)) 
   plt.ylabel('Cumulative variance (%)')
   plt.grid(True, alpha=0.3)
   plt.title('Cumulative variance captured by PCs')
   plt.savefig(outdir / 'cum_var.png')

   # Save 10 PCs to csv
   raw = np.loadtxt(outdir / 'proj_10.xvg', comments=['@','#','&'])  # shape: (2N, 2)
   n_frames = len(raw) // 10
   logger.debug(f"Loaded {n_frames} frames for PCs")
   pcs = []
   for i in range(10):
      pcs.append(raw[i*n_frames:(i+1)*n_frames, 1])  # PC i+1 is in column 1 of the corresponding block
   # Time from column 0 (same for both)
   time = raw[:n_frames, 0].flatten()  

   # Save to csv
   df_pc10 = pd.DataFrame({f'PC{i+1}': pcs[i] for i in range(10)})
   df_pc10.insert(0, 'Time(ps)', time)
   df_pc10.to_csv(outdir / 'proj_10.csv', index=False)

   # Plot PC1 vs PC2 colored by time
   plt.figure(figsize=(10, 8))
   scatter = plt.scatter(pcs[0], pcs[1], c=time, cmap='viridis', alpha=0.8, s=20)
   plt.colorbar(scatter, label='Time (ps)')
   plt.xlabel('PC1'); plt.ylabel('PC2')
   plt.title('PC1 vs PC2 colored by time')
   plt.grid(True, alpha=0.3)
   plt.tight_layout()
   plt.savefig(outdir / 'pc1_pc2.png', dpi=300)
   plt.close()

   # Pairwise scatter plots of first PCs
   fig, axes = plt.subplots(4, 4, figsize=(17, 15))
   pcs_top = ['PC1', 'PC2', 'PC3', 'PC4']

   for i, pc_x in enumerate(pcs_top):
      for j, pc_y in enumerate(pcs_top):
         ax = axes[i, j]
         
         if i == j: # Diagonal
               ax.scatter(df_pc10['Time(ps)'], df_pc10[pcs_top[i]], s=10, alpha=0.7, c='blue')
               ax.set_xlabel('Time (ps)')
               ax.set_ylabel(pcs_top[i])
               ax.set_title(f'{pcs_top[i]} vs Time')
         else:
               # Off-diagonal: scatter colored by Time
               scatter = ax.scatter(df_pc10[pc_y], df_pc10[pc_x], c=df_pc10['Time(ps)'], 
                                 cmap='viridis', s=10, alpha=0.7)
               ax.set_title(f'{pc_x} vs {pc_y}')
               ax.set_xlabel(pc_y)
               ax.set_ylabel(pc_x)

   plt.subplots_adjust(hspace=0.4, wspace=0.45)
   plt.colorbar(scatter, ax=axes.ravel().tolist(), label='Time (ps)', location='right', shrink=0.3, pad=0.03)
   plt.suptitle('Pairwise scatter plots of PCs colored by time', fontsize=16, y=0.94, x=0.45)
   plt.savefig(outdir / 'pc_pairgrid.png', dpi=300, bbox_inches='tight')
   plt.close()

   # RMSF plots of first 4 PCs
   rmsf = np.loadtxt(outdir / 'eigrmsf.xvg', comments=['@','#', '&'], usecols=[1])
   natoms = len(rmsf) // 10  # For 10 PCs
   fig, axes = plt.subplots(4, 1, figsize=(12, 15))
   for i, pc_x in enumerate(pcs_top):
      ax = axes[i]
      ax.plot(range(1, natoms + 1), rmsf[i*natoms:(i+1)*natoms], label=f'{pc_x} RMSF')
      ax.set_xticks(np.arange(0, natoms + 1, round(natoms/30)))
      ax.tick_params(axis='x', labelrotation=-90)
      ax.set_xlabel('CA atom')
      ax.set_ylabel('RMSF (nm)')
      ax.set_title(f'RMSF of {pc_x}', fontsize=14)
      ax.grid(True, alpha=0.3)
   plt.subplots_adjust(hspace=0.4)
   plt.savefig(outdir / 'pc_rmsf.png', bbox_inches='tight')
   plt.close()


if __name__ == "__main__":
    main()