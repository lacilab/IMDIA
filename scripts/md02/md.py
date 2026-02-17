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

def get_scaling(tpr_file, stage, total_cores = None):
    if total_cores is None:
        total_cores = get_total_cores()
    logger.debug(f"Determining scaling for {stage} stage with total available cores: {total_cores}")
    # Parse tpr with gmx dump
    dump = subprocess.check_output(['gmx', 'dump', '-s', tpr_file, '-nonr', '-quiet'], text=True)

    with open('dumptest.log', 'w') as log:
        log.write(f"Parsed tpr dump for {stage} stage:\n{dump}\n")
    
    # Parse separate fourier variables (your format)
    fourierspacing_match = re.search(r'fourierspacing\s*=\s*([\d.]+)', dump)
    fourier_nx_match = re.search(r'fourier-nx\s*=\s*(\d+)', dump)
    fourier_ny_match = re.search(r'fourier-ny\s*=\s*(\d+)', dump)
    fourier_nz_match = re.search(r'fourier-nz\s*=\s*(\d+)', dump)
    
    if all([fourierspacing_match, fourier_nx_match, fourier_ny_match, fourier_nz_match]):
        spacing = float(fourierspacing_match.group(1))
        gridx, gridy, gridz = map(int, [fourier_nx_match.group(1), 
                                       fourier_ny_match.group(1), 
                                       fourier_nz_match.group(1)])
        box_nm = min(gridx, gridy, gridz) * spacing
        logger.debug(f"Parsed fourier grid: {gridx}x{gridy}x{gridz} with spacing {spacing} nm, resulting in box size {box_nm} nm.")
    else:
        # Fallback: old compact format
        grid_match = re.search(r'fourier grid of (\d+)x(\d+)x(\d+).*spacing ([\d.]+)', dump)
        if grid_match:
            gridx, gridy, gridz, spacing = map(float, grid_match.groups())
            box_nm = min(gridx, gridy, gridz) * spacing
        else:
            box_nm = 10.0  # Default
        logger.debug(f"Parsed fourier grid (fallback): {gridx}x{gridy}x{gridz} with spacing {spacing} nm, resulting in box size {box_nm} nm.")
    
    # natoms 
    natoms_match = re.search(r'natoms\s*=\s*(\d+)', dump)
    natoms = int(natoms_match.group(1)) if natoms_match else 10000
    logger.debug(f"Parsed number of atoms: {natoms}")
    
    if stage == 'em':
        # Safe defaults based on box size
        if box_nm < 12:  
            ntmpi = 1      # Small box
        elif box_nm < 16:
            ntmpi = 2      # Medium box
        elif box_nm < 20:
            ntmpi = 4      # Medium box
        else:
            ntmpi = 8      # Large box
            
        if natoms < 50000:
            ntomp = min(total_cores // ntmpi // 4, 16)
        elif natoms < 200000:
            ntomp = min(total_cores // ntmpi // 2, 32)
        else:
            ntomp = min(total_cores // ntmpi, 64)
        return max(1, ntmpi), max(1, ntomp)
    
    elif stage in ['nvt', 'npt']:
        ntmpi = min(32, total_cores // 8)   
        ntomp = total_cores // ntmpi // 2 
        return max(1, ntmpi), max(1, ntomp)
    
    return 0, 0

def get_total_cores():
    if 'SLURM_CPUS_PER_TASK' in os.environ:
        return int(os.environ['SLURM_CPUS_PER_TASK'])
    elif 'SLURM_NTASKS' in os.environ:
        cpus_per_task = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
        return int(os.environ['SLURM_NTASKS']) * cpus_per_task
    
    # Fallback: physical cores
    return multiprocessing.cpu_count() // 2

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
        logger.debug(f"No custom parameter file provided for {filename}. Using default.")

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
    subprocess.run(cmd, check=True, text=True)
    # log_process_output(process, "gmx grompp EM")

    # Calculate resource use
    em_ntmpi, em_ntomp = get_scaling(f"{outdir}/em.tpr", 'em')
    logger.debug(f"Running EM with: gmx mdrun -ntmpi {em_ntmpi} -ntomp {em_ntomp}")

    cmd = ["gmx", "mdrun", "-v", "-deffnm", f"{outdir}/em", "-ntmpi", str(em_ntmpi), "-ntomp", str(em_ntomp), "-quiet"] + flag_nb_gpu
    subprocess.run(cmd, check=True, text=True)
    # log_process_output(process, "gmx mdrun EM")
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
    subprocess.run(cmd, check=True, text=True)
    # log_process_output(process, "gmx grompp NVT")
    nvt_ntmpi, nvt_ntomp = get_scaling(f"{outdir}/nvt.tpr", 'nvt')
    subprocess.run(["gmx", "mdrun", "-v", "-deffnm", f"{outdir}/nvt", "-ntmpi", str(nvt_ntmpi), "-ntomp", str(nvt_ntomp), "-quiet"] 
                   + flag_nb_gpu, check=True, text=True)
    # log_process_output(process, "gmx mdrun NVT")
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
    subprocess.run(cmd, check=True, text=True)
    # log_process_output(process, "gmx grompp NPT")
    npt_ntmpi, npt_ntomp = get_scaling(f"{outdir}/npt.tpr", 'npt')
    subprocess.run(["gmx", "mdrun", "-v", "-deffnm", f"{outdir}/npt", "-ntmpi", str(npt_ntmpi), "-ntomp", str(npt_ntomp), "-quiet"] 
                   + flag_nb_gpu, check=True, text=True)
    # log_process_output(process, "gmx mdrun NPT")
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
        subprocess.run(cmd, check=True, text=True)
        # log_process_output(process, "gmx mdrun MD (from checkpoint)")
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
        subprocess.run(cmd, check=True, text=True)
        # log_process_output(process, "gmx grompp MD")
        cmd = ["gmx", "mdrun", "-ntmpi", str(args.ntmpi), "-ntomp", str(args.ntomp)
               ] + flag_nb_gpu + flag_pme_gpu + flag_gpu_tasks + ["-v", "-deffnm", f"{outdir}/md", "-quiet"]
        logger.debug(f"Production MD command: {cmd}")
        subprocess.run(cmd, check=True, text=True)
        # log_process_output(process, "gmx mdrun MD")
    check_step("md", "Production MD", outdir)

# Create checkpoint file
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
checkpoint = os.path.join(outdir, f"check_02_md_{timestamp}.done")
with open(checkpoint, 'a'):
    os.utime(checkpoint, None)