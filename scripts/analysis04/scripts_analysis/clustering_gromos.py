#!/usr/bin/env python3

from pathlib import Path
import subprocess
from scripts.utils.helpers import setup_logger

# outdir OUTDIR   Output directory for GROMACS files   (required)
# tpr TPRFILE     TPR file for the system               (required)
# xtc XTCFILE     XTC file for the system               (required)

def main(logfile: str, outdir: Path, tpr: Path, xtc: Path):
    # Initialize logger
    logger = setup_logger("analysis", logfile) 

    # If outdir does not exist, create it
    if not outdir.is_dir():
        outdir.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Created output directory: {outdir}")
    else:
        logger.debug(f"Output directory already exists: {outdir}")

    # Check if tpr and xtc files exist
    if not tpr.is_file():
        logger.error(f"TPR file not found: {tpr}")
        return
    else:
        logger.debug(f"Found TPR file: {tpr}")

    if not xtc.is_file():
        logger.error(f"XTC file not found: {xtc}")
        return
    else:
        logger.debug(f"Found XTC file: {xtc}")

    # Run gromos clustering
    cmd = ["gmx", "cluster", "-f", str(xtc), "-s", str(tpr), "-method", "gromos", "-cutoff", "0.1",
            "-o", str(outdir / "gromos_clusters.xpm"), "-g", str(outdir / "gromos_clusters.log"),
            "-clndx", str(outdir / "gromos_clusters.ndx"), "-dist", str(outdir / "gromos_clusters_rmsd_dist.xvg"),
            "-clid", str(outdir / "gromos_clust_ids.xvg"), "-sz", str(outdir / "gromos_clust_size.xvg"), "-quiet"]
    subprocess.run(cmd, check=True, text=True, input="C-alpha\nC-alpha\n")
    if (outdir / "gromos_clusters.xpm").is_file():
        logger.info(f"GROMOS clustering completed successfully. Output files saved in {outdir}")

    # Create EPS file from XPM
    cmd = ["gmx", "xpm2ps", "-f", str(outdir / "gromos_clusters.xpm"), "-o", str(outdir / "gromos_clusters.eps"), "-quiet"]
    subprocess.run(cmd, check=True, text=True)

    # Convert EPS to PNG
    cmd = ["gs", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r300", "-dEPSCrop",
           "-sOutputFile=" + str(outdir / "gromos_clusters.png"),
           str(outdir / "gromos_clusters.eps")]
    subprocess.run(cmd, check=True, text=True)
    if (outdir / "gromos_clusters.png").is_file():
        logger.info(f"GROMOS cluster matrix plot saved as PNG: {outdir / 'gromos_clusters.png'}")



if __name__ == "__main__":
    main()