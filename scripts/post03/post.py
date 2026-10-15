from datetime import datetime
import os
from pathlib import Path
import subprocess
import sys
from scripts.utils.helpers import setup_logger
from scripts.utils.group_index import group_index

def main(logfile: Path, repdir: Path, dna: bool = False, rna: bool = False, phos: bool = False):
    """
    Post-process MD results to make them ready for analysis and visualization.

    Steps:
       1.Remove PBC artifacts.
       2.Cluster complex.
       3.Fit the system to reference structure (using Backbone).
       4.Cleanup intermediate files.
       5.Output post-processed trajectories and PDBs at frame 10.

    Arguments:
       logfile : Path to the log file.
       repdir  : Directory containing the MD results for this replicate.
       dna     : True if the system contains DNA.
       rna     : True if the system contains RNA.
       phos    : True if the system contains phosphorylated residues.
    """

    # Initialize logger
    logger = setup_logger("post", logfile)

    def log_process_output(process, step_name):
        if process.returncode != 0:
            logger.error(f"{step_name} error:\n {process.stderr}")
        elif process.stdout:
            logger.debug(f"{step_name} output:\n {process.stdout}")  

    checkpoint_files = list(repdir.glob("check_02_post_*.done"))
    if len(checkpoint_files) > 0:
        logger.warning("Checkpoint file found. Skipping post-processing.")
        logger.warning(f"To rerun these steps, please remove the checkpoint file(s): {repdir}/check_02_post_*.done")
        sys.exit(0)

    # Get main group name
    group = group_index(gmxdir=repdir, logfile=logfile, logname='post', structurefile=repdir / "md.gro", dna=dna, rna=rna, phos=phos)

    ##### POST-PROCESSING STEPS ######
    # Remove PBC artifacts
    cmd = ["gmx", "trjconv", 
           "-s", repdir / "md.tpr", 
           "-f", repdir / "md.xtc", 
           "-o", repdir / "md_noPBC_mol_com.xtc", 
           "-pbc", "mol", 
           "-ur", "compact", 
           "-quiet", "-n", repdir / "index.ndx"]
    process = subprocess.run(cmd, input="System\n", text=True, check=True, capture_output=True)
    log_process_output(process, "gmx trjconv - remove PBC")

    # Cluster to avoid monomers jumping 
    cmd = ["gmx", "trjconv", 
           "-s", repdir / "md.tpr", 
           "-f", repdir / "md_noPBC_mol_com.xtc", 
           "-o", repdir / "md_noPBC_mol_com_cluster.xtc", 
           "-pbc", "cluster", 
           "-quiet", "-n", repdir / "index.ndx"]
    input_str = f"{group}\nSystem\n"
    process = subprocess.run(cmd, input=input_str, text=True, check=True, capture_output=True)
    log_process_output(process, "gmx trjconv - cluster")

    # Fit the system to reference structure (using Backbone)
    cmd = ["gmx", "trjconv",
           "-s", repdir / "md.tpr", 
           "-f", repdir / "md_noPBC_mol_com_cluster.xtc", 
           "-o", repdir / "md_noPBC_mol_com_cluster_fit.xtc", 
           "-fit", "rot+trans", 
           "-quiet", 
           "-n", repdir / "index.ndx"]
    input_str = f"Backbone\nSystem\n"
    process = subprocess.run(cmd, input=input_str, text=True, check=True, capture_output=True)
    log_process_output(process, "gmx trjconv - fit")

    # Clean up intermediate files
    os.remove(repdir / "md_noPBC_mol_com.xtc")
    os.remove(repdir / "md_noPBC_mol_com_cluster.xtc")

    # Rename post-processed trajectory
    os.rename(repdir / "md_noPBC_mol_com_cluster_fit.xtc", repdir / "md_post.xtc")

    # Output .pdb files at frame 10
    cmd = ["gmx", "trjconv", 
           "-s", repdir / "md.tpr", 
           "-f", repdir / "md_post.xtc", 
           "-o", repdir / "md_sys.pdb", 
           #"-pbc", "mol", 
           #"-ur", "compact", 
           "-dump", "10", 
           "-quiet", "-n", repdir / "index.ndx"]
    process = subprocess.run(cmd, input="System\n", check=True, text=True, capture_output=True)
    log_process_output(process, "gmx trjconv - output pdb")

    cmd = ["gmx", "trjconv", 
           "-s", repdir / "md.tpr", 
           "-f", repdir / "md_post.xtc", 
           "-o", repdir / "md_prot.pdb", 
           #"-pbc", "mol", 
           #"-ur", "compact", 
           "-dump", "10", 
           "-quiet", "-n", repdir / "index.ndx"]
    process = subprocess.run(cmd, input=f"{group}\n", check=True, text=True, capture_output=True)
    log_process_output(process, "gmx trjconv - output pdb protein complex")

    # Create checkpoint file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    checkpoint = repdir / f"check_03_post_{timestamp}.done"
    checkpoint.parent.mkdir(exist_ok=True, parents=True)
    with checkpoint.open('w') as f:
        f.write(f"Post-processing completed at {timestamp}\n")
        logger.info(f"Post-processing completed successfully. Checkpoint created: {checkpoint}")
        return


if __name__ == "__main__":
    main()