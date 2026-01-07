#!/usr/bin/env python
# -*- coding: utf-8 -*-

import datetime
import sys
import os
import glob
import re
import multiprocessing
import subprocess
import shutil
import argparse
from scripts.utils.helpers import setup_logger
from scripts.utils.group_index import group_index

parser = argparse.ArgumentParser(description="Perform MD simulation steps: minimization, NVT, NPT, and production MD.")
parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
parser.add_argument('--gmxdir', '-g', dest='gmxdir', default='.', help='Directory containing GROMACS input files (ions.gro, topol.top, etc.)')
parser.add_argument('--repn', '-n', dest='repn', type=int, default=None, help='Replicate number (used to create output directory)')
parser.add_argument('--dna', dest='with_dna', action='store_true', help='Use if the system contains DNA. Adds DNA.bsc1 force field.')
parser.add_argument('--phos', dest='with_phosphate', action='store_true', help='Use if the system contains phosphorylated residues.')
parser.add_argument('--ntmpi', dest='ntmpi', type=int, default=0, help='Number of MPI ranks (default: auto-detected)')
parser.add_argument('--ntomp', dest='ntomp', type=int, default=0, help='Number of OpenMP threads per rank (default: auto-detected)')
parser.add_argument('--nb', dest='nb', type=str, default='gpu', help='Calculate non-bonded interactions on: auto, cpu, gpu')
parser.add_argument('--pme', dest='pme', type=str, default='gpu', help='Perform PME calculations on: auto, cpu, gpu')
parser.add_argument('--npme', dest='npme', type=int, default=-1, help='Number of separate ranks to be used for PME, -1 is guess (let GROMACS decide)')
parser.add_argument('--ngpu', dest='ngpu', type=int, default=None, help='Number of GPUs to use (default: auto-detected)')
parser.add_argument('--gputasks', dest='gputasks', type=str, default=None, help='Sequence of GPU IDs to use i.e. 0011 (default: auto-detected)')
parser.add_argument('--em', dest='em', type=str, default=None, help='Energy minimization parameter file (default: minim.mdp)')
parser.add_argument('--nvt', dest='nvt', type=str, default=None, help='NVT equilibration parameter file (default: nvt.mdp)')    
parser.add_argument('--npt', dest='npt', type=str, default=None, help='NPT equilibration parameter file (default: npt.mdp)')
parser.add_argument('--md', dest='md', type=str, default=None, help='Production MD parameter file (default: md.mdp)')

# print help if no arguments are provided
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(0)

args = parser.parse_args()

##### INITIALIZE LOG
logger = setup_logger("md", args.logfile)

##### HELPER FUNCTIONS #####

# helper to set parameter files. args: <custom_path_var> <filename>
# if -em, --nvt, --npt, --md options are provided, copy those files instead of the default ones to $outdir. Else, copy the default files
def param_mdp(filename, outdir, custom_path = None):
    if custom_path is not None:
        if os.path.isfile(custom_path):
            shutil.copy(custom_path, f"{outdir}/{filename}")
            logger.info(f"Using custom parameter file {custom_path} for {filename}.")
        else:
            logger.warning(f"Custom parameter file {custom_path} not found. Using default {filename}.")
            shutil.copy(f"scripts/md02/{filename}", f"{outdir}/{filename}")
    else:
        shutil.copy(f"scripts/md02/{filename}", f"{outdir}/{filename}")

# Helper to read last N lines of the file (adapted from GeeksforGeeks)
def LastNlines(fname, N):
    lines = []
    with open(fname) as file:
        for line in (file.readlines()[-N:]):
            lines.append(line)
    return lines

def check_finished(step_short, outdir):
    if os.path.isfile(f"{outdir}/{step_short}.gro") and os.path.isfile(f"{outdir}/{step_short}.log") and "Finished mdrun" in " ".join(LastNlines(f"{outdir}/{step_short}.log", 20)):
        return True
    else:
        return False
    
def log_process_output(process, step_name):
    if process.returncode != 0:
        logger.error(f"{step_name} error:\n {process.stderr}")
    elif process.stdout:
        logger.debug(f"{step_name} output:\n {process.stdout}")   


# helper to check the step outputs.
def check_step(step_short, step_name, outdir):
    if check_finished(step_short, outdir):
        logger.info(f"{step_name} completed successfully.")
    else:
        logger.error(f"{step_name} failed. Check {outdir}/{step_short}.log for details.")
        sys.exit(1)


##### PARAMETER CHECKS #####

logger.debug(f"Received positional argument gmxdir: {args.gmxdir} and replicate number: {args.repn}")

if not os.path.isdir(args.gmxdir):
    logger.error(f"GROMACS directory {args.gmxdir} does not exist.")
    sys.exit(1)

parent_dir = os.path.dirname(args.gmxdir)

if args.repn == None:
    # List all subdirectories in parent matching rep*
    rep_dirs = [
        d for d in os.listdir(parent_dir)
        if os.path.isdir(os.path.join(parent_dir, d)) and re.match(r'rep\d+', d)
    ]
    # Extract numbers
    rep_numbers = [int(re.search(r'rep(\d+)', d).group(1)) for d in rep_dirs]
    # Sort numerically
    rep_numbers.sort()
    # Determine next replication number
    if not rep_numbers:
        n = 1
    else:
        n = rep_numbers[-1] + 1
    logger.debug("Next rep number:", n)
else:
    n = args.repn

outdir = f"{parent_dir}/rep{n}/gromacs"
# check if outdir exists, if not create it
if not os.path.isdir(outdir):
    os.makedirs(outdir)
    logger.info(f"Created output directory: {outdir}")

# Check for checkpoint files
if glob.glob(f"{outdir}/check_02_md_*.done"):
    logger.warning("Checkpoint file found. Skipping MD steps: minimization, NVT, NPT and production MD.")
    logger.debug(f"To rerun these steps, please remove the checkpoint file(s): {outdir}/check_02_md_*.done")
    sys.exit(0)

# TODO Not necessary, let GROMACS detect it
# Check number of GPUs if not provided
# if args.ngpu is None:
#     try:
#         # Check if 'nvidia-smi' command exists
#         subprocess.run(['nvidia-smi'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
#         # Count the number of GPUs
#         result = subprocess.run(['nvidia-smi', '-L'], stdout=subprocess.PIPE, text=True)
#         NGPUS = len(result.stdout.strip().splitlines())
#     except (subprocess.CalledProcessError, FileNotFoundError):
#         # nvidia-smi not found or failed
#         NGPUS = 0
# else:
#     NGPUS = args.ngpu

# TODO Not necessary, set default NTMPI = 0 and NTOMP = 0 in argparse and let GROMACS decide
# Auto-detect NTMPI and NTOMP if not provided
# SLURM_NTASKS = os.environ.get('SLURM_NTASKS')
# SLURM_CPUS_PER_TASK = os.environ.get('SLURM_CPUS_PER_TASK')
# 
# if SLURM_NTASKS and SLURM_CPUS_PER_TASK:
#     NTMPI = int(SLURM_NTASKS)
#     NTOMP = int(SLURM_CPUS_PER_TASK)
# else:
#     # fallback for local machine
#     NTOMP = multiprocessing.cpu_count()
#     NTMPI = 1
# 
# if args.ntmpi is not None:
#     NTMPI = args.ntmpi
# if args.ntomp is not None:
#     NTOMP = args.ntomp

##### RESOURCE SETUP #####
logger.debug(f"Running {args.ntmpi} MPI ranks x {args.ntomp} threads per rank.")

# Initialize empty flags
flag_nb_gpu = []
flag_pme_gpu = []
flag_gpu_tasks = []
GPUTASKS = ""

# Assign flags based on resources
if args.ngpu is not None:
    flag_nb_gpu = ["-nb", args.nb]
    flag_pme_gpu = ["-pme", args.pme, "-npme", str(args.npme)]
    if args.gputasks is not None:
        GPUTASKS = args.gputasks
        logger.debug(f"Using specified GPU list: {GPUTASKS}")
    else:
        logger.debug("No GPU task list specified. Generating GPU task list automatically.")
        for i in range(1, args.ntmpi+1):
            gpu = i % args.ngpu if args.ngpu > 0 else 0
            GPUTASKS += str(gpu)
    flag_gpu_tasks = ["-gputasks", GPUTASKS]
    logger.debug(f"Using {args.ngpu} GPU(s) with gputasks: {GPUTASKS}")
else:
    logger.debug("No GPUs specified. Detecting available GPUs and assigning the default behavior.")


##### PARAMETER FILES SETUP #####
# copy parameter files to outdir
# parameter files
param_mdp('minim.mdp', outdir, args.em)
param_mdp('nvt.mdp', outdir, args.nvt)
param_mdp('npt.mdp', outdir, args.npt)
param_mdp('md.mdp', outdir, args.md)

# check outdir exists and contains minim.mdp, nvt.mdp, npt.mdp, md.mdp
if not os.path.isdir(outdir):
    logger.error(f"Output directory {outdir} does not exist.")
    sys.exit(1)

for mdp_file in ['minim.mdp', 'nvt.mdp', 'npt.mdp', 'md.mdp']:
    if not os.path.isfile(f"{outdir}/{mdp_file}"):
        logger.error(f"Parameter file {mdp_file} not found in output directory {outdir}.")
        sys.exit(1)

logger.debug(f"Working directory: {os.getcwd()}, gmxdir: {args.gmxdir}, outdir: {outdir}")

##### MD STEPS #####

### 1. Energy Minimization
if check_finished("em", outdir):
    logger.warning("em.log indicates minimization already completed and em.gro exists. Skipping energy minimization.")
else:
    logger.debug("Starting energy minimization...")
    cmd = ["gmx", "grompp",
           "-f", f"{outdir}/minim.mdp",
           "-c", f"{args.gmxdir}/ions.gro",
           "-r", f"{args.gmxdir}/ions.gro",
           "-p", f"{args.gmxdir}/topol.top",
           "-o", f"{outdir}/em.tpr",
           "-po", f"{outdir}/em_out.mdp",
           "-quiet", "-maxwarn", "1"]
    process = subprocess.run(cmd, check=True, capture_output=True, text=True)
    log_process_output(process, "gmx grompp EM")

    cmd = ["gmx", "mdrun", "-v", "-deffnm", f"{outdir}/em", "-quiet"] + flag_nb_gpu
    process = subprocess.run(cmd, check=True, capture_output=True, text=True)
    log_process_output(process, "gmx mdrun EM")
    check_step("em", "Energy Minimization", outdir)

### 2. Create groups
group_index(gmxdir=outdir, logfile=args.logfile, logname='md', structurefile=f"{outdir}/em.gro", dna=args.with_dna, phos=args.with_phosphate)

### 3. NVT Equilibration
if check_finished("nvt", outdir):
    logger.warning("nvt.log indicates NVT equilibration already completed and nvt.gro exists. Skipping NVT equilibration.")
else:
    logger.debug("Starting NVT equilibration...")
    cmd = ["gmx", "grompp",
            "-f", f"{outdir}/nvt.mdp",
            "-c", f"{outdir}/em.gro",
            "-r", f"{outdir}/em.gro",
            "-p", f"{args.gmxdir}/topol.top",
            "-o", f"{outdir}/nvt.tpr",
            "-po", f"{outdir}/nvt_out.mdp",
            "-n", f"{outdir}/index.ndx", 
            "-quiet", "-maxwarn", "1"]
    process = subprocess.run(cmd, check=True, text=True, capture_output=True)
    log_process_output(process, "gmx grompp NVT")
    process = subprocess.run(["gmx", "mdrun", "-v", "-deffnm", f"{outdir}/nvt", "-quiet"] + flag_nb_gpu)
    log_process_output(process, "gmx mdrun NVT")
    check_step("nvt", "NVT Equilibration", outdir)

### 4. NPT Equilibration
if check_finished("npt", outdir):
    logger.warning("npt.log indicates NPT equilibration already completed and npt.gro exists. Skipping NPT equilibration.")
else:
    logger.debug("Starting NPT equilibration...")
    cmd = ["gmx", "grompp",
           "-f", f"{outdir}/npt.mdp",
           "-c", f"{outdir}/nvt.gro",
           "-r", f"{outdir}/nvt.gro",
           "-t", f"{outdir}/nvt.cpt",
           "-p", f"{args.gmxdir}/topol.top",
           "-o", f"{outdir}/npt.tpr",
           "-po", f"{outdir}/npt_out.mdp",
           "-n", f"{outdir}/index.ndx",
           "-quiet", "-maxwarn", "1"]
    process = subprocess.run(cmd, check=True, text=True, capture_output=True)
    log_process_output(process, "gmx grompp NPT")
    process = subprocess.run(["gmx", "mdrun", "-v", "-deffnm", f"{outdir}/npt", "-quiet"] + flag_nb_gpu)
    log_process_output(process, "gmx mdrun NPT")
    check_step("npt", "NPT Equilibration", outdir)

### 5. Production MD
if check_finished("md", outdir):
    logger.warning("md.log indicates production MD already completed and md.gro exists. Skipping production MD.")
else:
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    logfile = os.path.join(outdir, f"md_{timestamp}.log")
    if os.path.isfile(f"{outdir}/md.cpt"):
        logger.debug(f"Resuming production MD from checkpoint {outdir}/md.cpt...")
        cmd = ["gmx", "mdrun", 
            "-ntmpi", str(args.ntmpi), 
            "-ntomp", str(args.ntomp)] + flag_nb_gpu + flag_pme_gpu + flag_gpu_tasks + [
                "-cpi", f"{outdir}/md.cpt", "-append", "-v", "-deffnm", f"{outdir}/md", "-quiet"]
        process = subprocess.run(cmd, capture_output=True, text=True, check=True)
        log_process_output(process, "gmx mdrun MD (from checkpoint)")
    else:
        logger.debug("Starting production MD...")
        cmd = ["gmx", "grompp",
               "-f", f"{outdir}/md.mdp",
               "-c", f"{outdir}/npt.gro",
               "-r", f"{outdir}/npt.gro",
               "-t", f"{outdir}/npt.cpt",
               "-p", f"{args.gmxdir}/topol.top",
               "-o", f"{outdir}/md.tpr",
                "-po", f"{outdir}/md_out.mdp",
               "-n", f"{outdir}/index.ndx",
               "-quiet", "-maxwarn", "1"]
        process = subprocess.run(cmd, check=True, text=True, capture_output=True)
        log_process_output(process, "gmx grompp MD")
        cmd = ["gmx", "mdrun", "-ntmpi", str(args.ntmpi), "-ntomp", str(args.ntomp)
               ] + flag_nb_gpu + flag_pme_gpu + flag_gpu_tasks + ["-v", "-deffnm", f"{outdir}/md", "-quiet"]
        logger.debug(f"Production MD command: {cmd}")
        process = subprocess.run(cmd, capture_output=True, text=True, check=True)
        log_process_output(process, "gmx mdrun MD")
    check_step("md", "Production MD", outdir)

# Create checkpoint file
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
checkpoint = os.path.join(outdir, f"check_02_md_{timestamp}.done")
with open(checkpoint, 'a'):
    os.utime(checkpoint, None)