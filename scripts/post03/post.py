
import argparse
from datetime import datetime
import glob
import os
import shutil
import subprocess
import sys
from scripts.utils.helpers import setup_logger



def main():
    parser = argparse.ArgumentParser(description="Perform MD simulation steps: minimization, NVT, NPT, and production MD.")
    parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
    parser.add_argument('--gmxdir', '-g', dest='gmxdir', default='.', help='Directory containing GROMACS input files (ions.gro, topol.top, etc.)')
    parser.add_argument('--dna', dest='with_dna', action='store_true', help='Use if the system contains DNA. Adds DNA.bsc1 force field.')
    parser.add_argument('--phos', dest='with_phosphate', action='store_true', help='Use if the system contains phosphorylated residues.')

    # print help if required arguments are not provided
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Initialize logger
    logger = setup_logger("post", args.logfile)

    if glob.glob(f"{args.gmxdir}/check_02_post_*.done"):
        logger.warning("Checkpoint file found. Skipping post-processing.")
        print(f"To rerun these steps, please remove the checkpoint file(s): {args.gmxdir}/check_02_post_*.done")
        sys.exit(0)

    # Construct group name
    group = "Protein"
    if args.with_dna:
        group += "_DNA"
    if args.with_phosphate:
        group += "_P"

    def log_process_output(process, step_name):
        if process.returncode != 0:
            logger.error(f"{step_name} error:\n {process.stderr}")
        elif process.stdout:
            logger.debug(f"{step_name} output:\n {process.stdout}")  

    ##### POST-PROCESSING STEPS ######
    # center molecule and put all atoms at the closest distance from center of the box
    # TODO check manually with active and inactive structures, added -center from the previous version
    cmd = ["gmx", "trjconv", 
           "-s", "md.tpr", 
           "-f", "md.xtc", 
           "-o", "md_noPBC_mol_com.xtc", 
           "-pbc", "mol", 
           "-ur", "compact", 
           "-center", 
           "-quiet", "-n", "index.ndx"]
    process = subprocess.run(cmd, input="0\n", check=True, capture_output=True)
    log_process_output(process, "gmx trjconv - remove PBC")

    # Cluster to avoid monomers jumping 
    cmd = ["gmx", "trjconv", 
           "-s", "md.tpr", 
           "-f", "md_noPBC_mol_com.xtc", 
           "-o", "md_noPBC_mol_com_cluster.xtc", 
           "-pbc", "cluster", 
           "-quiet", "-n", "index.ndx"]
    input_str = f"{group}\n0\n"
    process = subprocess.run(cmd, input=input_str, text=True, check=True, capture_output=True)
    log_process_output(process, "gmx trjconv - cluster")

    # Fit the system to reference structure (using Backbone)
    cmd = ["gmx", "trjconv",
           "-s", "md.tpr", 
           "-f", "md_noPBC_mol_com_cluster.xtc", 
           "-o", "md_noPBC_mol_com_cluster_fit.xtc", 
           "-fit", "rot+trans", 
           "-quiet", 
           "-n", "index.ndx"]
    input_str = f"Backbone\n0\n"
    process = subprocess.run(cmd, input=input_str, text=True, check=True, capture_output=True)
    log_process_output(process, "gmx trjconv - fit")

    # Clean up intermediate files
    shutil.remove("md_noPBC_mol_com.xtc")
    shutil.remove("md_noPBC_mol_com_cluster.xtc")

    # Output .pdb files for original trajectory reference (t=0 of raw simulation)
    cmd = ["gmx", "trjconv", 
           "-s", "md.tpr", 
           "-f", "md.xtc", 
           "-o", "md_sys.pdb", 
           "-pbc", "mol", 
           "-ur", "compact", 
           "-dump", "0", 
           "-quiet", "-n", "index.ndx"]
    process = subprocess.run(cmd, input="0\n", check=True, capture_output=True)
    log_process_output(process, "gmx trjconv - output pdb")

    process = subprocess.run(cmd, input=f"{group}\n", check=True, capture_output=True)
    log_process_output(process, "gmx trjconv - output pdb protein complex")

    # Create checkpoint file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    checkpoint = os.path.join(args.gmxdir, f"check_03_post_{timestamp}.done")
    with open(checkpoint, 'w') as f:
        f.write(f"Post-processing completed at {timestamp}\n")

    logger.info("Post-processing completed successfully. Checkpoint created: " + checkpoint)


if __name__ == "__main__":
    main()