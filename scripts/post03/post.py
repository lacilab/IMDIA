
import argparse
from datetime import datetime
import glob
import os
import subprocess
import sys
from scripts.utils.helpers import setup_logger



def main():
    parser = argparse.ArgumentParser(description="Perform MD simulation steps: minimization, NVT, NPT, and production MD.")
    parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
    parser.add_argument('--repdir', '-n', dest='repdir', type=str, default=None, help='Directory containing the MD results for this replicate.')
    parser.add_argument('--dna', dest='with_dna', action='store_true', help='Use if the system contains DNA. Adds DNA.bsc1 force field.')
    parser.add_argument('--phos', dest='with_phosphate', action='store_true', help='Use if the system contains phosphorylated residues.')

    # print help if required arguments are not provided
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Initialize logger
    logger = setup_logger("post", args.logfile)

    if glob.glob(f"{args.repdir}/check_02_post_*.done"):
        logger.warning("Checkpoint file found. Skipping post-processing.")
        print(f"To rerun these steps, please remove the checkpoint file(s): {args.repdir}/check_02_post_*.done")
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
    # Remove PBC artifacts
    cmd = ["gmx", "trjconv", 
           "-s", f"{args.repdir}/md.tpr", 
           "-f", f"{args.repdir}/md.xtc", 
           "-o", f"{args.repdir}/md_noPBC_mol_com.xtc", 
           "-pbc", "mol", 
           "-ur", "compact", 
           "-quiet", "-n", f"{args.repdir}/index.ndx"]
    process = subprocess.run(cmd, input="System\n", text=True, check=True, capture_output=True)
    log_process_output(process, "gmx trjconv - remove PBC")

    # Cluster to avoid monomers jumping 
    cmd = ["gmx", "trjconv", 
           "-s", f"{args.repdir}/md.tpr", 
           "-f", f"{args.repdir}/md_noPBC_mol_com.xtc", 
           "-o", f"{args.repdir}/md_noPBC_mol_com_cluster.xtc", 
           "-pbc", "cluster", 
           "-quiet", "-n", f"{args.repdir}/index.ndx"]
    input_str = f"{group}\nSystem\n"
    process = subprocess.run(cmd, input=input_str, text=True, check=True, capture_output=True)
    log_process_output(process, "gmx trjconv - cluster")

    # Fit the system to reference structure (using Backbone)
    cmd = ["gmx", "trjconv",
           "-s", f"{args.repdir}/md.tpr", 
           "-f", f"{args.repdir}/md_noPBC_mol_com_cluster.xtc", 
           "-o", f"{args.repdir}/md_noPBC_mol_com_cluster_fit.xtc", 
           "-fit", "rot+trans", 
           "-quiet", 
           "-n", f"{args.repdir}/index.ndx"]
    input_str = f"Backbone\nSystem\n"
    process = subprocess.run(cmd, input=input_str, text=True, check=True, capture_output=True)
    log_process_output(process, "gmx trjconv - fit")

    # Clean up intermediate files
    os.remove(f"{args.repdir}/md_noPBC_mol_com.xtc")
    os.remove(f"{args.repdir}/md_noPBC_mol_com_cluster.xtc")

    # Rename post-processed trajectory
    os.rename(f"{args.repdir}/md_noPBC_mol_com_cluster_fit.xtc", f"{args.repdir}/md_post.xtc")

    # Output .pdb files at frame 10
    cmd = ["gmx", "trjconv", 
           "-s", f"{args.repdir}/md.tpr", 
           "-f", f"{args.repdir}/md_post.xtc", 
           "-o", f"{args.repdir}/md_sys.pdb", 
           #"-pbc", "mol", 
           #"-ur", "compact", 
           "-dump", "10", 
           "-quiet", "-n", f"{args.repdir}/index.ndx"]
    process = subprocess.run(cmd, input="System\n", check=True, text=True, capture_output=True)
    log_process_output(process, "gmx trjconv - output pdb")

    cmd = ["gmx", "trjconv", 
           "-s", f"{args.repdir}/md.tpr", 
           "-f", f"{args.repdir}/md_post.xtc", 
           "-o", f"{args.repdir}/md_prot.pdb", 
           #"-pbc", "mol", 
           #"-ur", "compact", 
           "-dump", "10", 
           "-quiet", "-n", f"{args.repdir}/index.ndx"]
    process = subprocess.run(cmd, input=f"{group}\n", check=True, text=True, capture_output=True)
    log_process_output(process, "gmx trjconv - output pdb protein complex")

    # Create checkpoint file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    checkpoint = os.path.join(args.repdir, f"check_03_post_{timestamp}.done")
    with open(checkpoint, 'w') as f:
        f.write(f"Post-processing completed at {timestamp}\n")

    logger.info("Post-processing completed successfully. Checkpoint created: " + checkpoint)


if __name__ == "__main__":
    main()