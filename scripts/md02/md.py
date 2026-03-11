#!/usr/bin/env python
# -*- coding: utf-8 -*-

import datetime
import sys
import os
import re
import multiprocessing
import subprocess
import shutil
from sympy import isprime
import numpy as np
from pathlib import Path
import math
import matplotlib.pyplot as plt
from scripts.utils.helpers import setup_logger
from scripts.utils.group_index import group_index

##### HELPER FUNCTIONS #####

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
def param_mdp(logger, filename: str, outdir: Path, custom_path: Path = None):
    if custom_path is not None:
        if custom_path.is_file():
            shutil.copy(custom_path, outdir / filename)
            logger.info(f"Using custom parameter file {custom_path} for {filename}.")
        else:
            logger.warning(f"Custom parameter file {custom_path} not found. Using default {filename}.")
            shutil.copy(Path(f"scripts/md02/{filename}"), outdir / filename)
    else:
        shutil.copy(Path(f"scripts/md02/{filename}"), outdir / filename)
        logger.debug(f"No custom parameter file provided for {filename}. Using default.")

# Helper to read last N lines of the file (adapted from GeeksforGeeks)
def LastNlines(fname, N):
    lines = []
    with open(fname) as file:
        for line in (file.readlines()[-N:]):
            lines.append(line)
    return lines

def check_finished(step_short: str, outdir: Path):
    if (outdir / f"{step_short}.gro").is_file() and (outdir / f"{step_short}.log").is_file() and "Finished mdrun" in " ".join(LastNlines(outdir / f"{step_short}.log", 20)):
        return True
    else:
        return False
    
def log_process_output(process, step_name, logger):
    if process.returncode != 0:
        logger.error(f"{step_name} error:\n {process.stderr}")
    elif process.stdout:
        logger.debug(f"{step_name} output:\n {process.stdout}")   

# to check the step outputs.
def check_step(logger, step_short: str, step_name: str, outdir: Path):
    if check_finished(step_short, outdir):
        logger.info(f"{step_name} completed.")
    else:
        logger.error(f"{step_name} failed. Check {outdir}/{step_short}.log for details.")
        sys.exit(1)

def get_grid_size(logger, stage:str, dump):    
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
    return box_nm

def get_natoms(logger, dump):
    natoms_match = re.search(r'natoms\s*=\s*(\d+)', dump)
    natoms = int(natoms_match.group(1)) if natoms_match else 10000
    logger.debug(f"Parsed number of atoms: {natoms}")
    return natoms

def get_mincellsize(logger, logfile: Path):
    """Extract minimum cell size from EM log"""

    if not logfile.is_file():
        logger.warning(f"EM log file '{logfile}' not found. Cannot extract minimum cell size for scaling decisions.")
        return None

    with logfile.open('r') as f:
        content = f.read()
    
    # Look for the DD info line (appears even with ntmpi=1)
    min_bonded = re.search(r'Minimum cell size due to bonded interactions: ([\d.]+)', content)
    min_disp = re.search(r'Minimum cell size due to atom displacement: ([\d.]+)', content)
    # Safe extraction with fallbacks
    bonded_size = float(min_bonded.group(1)) if min_bonded else None
    disp_size = float(min_disp.group(1)) if min_disp else None
    # Take max of available values
    available_sizes = [s for s in [bonded_size, disp_size] if s is not None]
    if available_sizes:
        mincellsize = max(available_sizes)
        logger.debug(f"Extracted minimum cell size from log: {mincellsize} nm")
        return mincellsize
    
    # Fallback patterns
    fallback_matches = [
        r'cell size limit.*?([\d.]+)',
        r'domain decomposition.*?([\d.]+)'
    ]
    for pattern in fallback_matches:
        match = re.search(pattern, content)
        if match:
            logger.debug(f"Extracted minimum cell size from fallback pattern: {match.group(1)} nm")
            return float(match.group(1))
    
    return None

def get_scaling(logger, tpr_file: Path, stage: str, total_cores = None):
    if total_cores is None:
        total_cores = get_total_cores()
    logger.debug(f"Determining scaling for {stage} stage with total available cores: {total_cores}")

    # Parse tpr with gmx dump
    dump = subprocess.check_output(['gmx', 'dump', '-s', tpr_file, '-nonr', '-quiet'], text=True)
    box_nm = get_grid_size(logger, stage, dump)
    natoms = get_natoms(logger, dump)
    mincellsize = get_mincellsize(logger, tpr_file.parent / "em.log")

    # For EM, run with conservative defaults to avoid decomposition errors
    if stage == "em" or mincellsize is None:
        logger.debug(f"Minimum cell size not found. Using safe defaults for scaling based on box size and number of atoms for step {stage}.")
        # Safe defaults based on box size
        if box_nm < 16:  
            ntmpi = 1      # Small box
        elif box_nm < 20:
            ntmpi = 2      # Medium box
        elif box_nm < 24:
            ntmpi = 4      # Medium box
        else:
            ntmpi = 8      # Large box
        ntmpi = min(total_cores, ntmpi)    
        if natoms < 50000:
            ntomp = min(total_cores // ntmpi // 4, 16)
        elif natoms < 200000:
            ntomp = min(total_cores // ntmpi // 2, 32)
        else:
            ntomp = min(total_cores // ntmpi, 64)
        return max(1, ntmpi), max(1, ntomp)
    # If minimum cell size is available, use it to calculate safe scaling with better resource use
    else:
       logger.debug(f"Minimum cell size found: {mincellsize} nm.")
       ntmpi = min(max(1, int(box_nm // mincellsize)), total_cores, 16)
       if ntmpi > 7 and isprime(ntmpi):
        ntmpi -= 1
       ntomp = max(1, total_cores // ntmpi)
       logger.debug(f"Calculated ntmpi {ntmpi} and ntomp {ntomp} based on box size {box_nm} and min cell size {mincellsize} for step {stage}.")
       return ntmpi, ntomp

# parser = argparse.ArgumentParser(description="Perform MD simulation steps: minimization, NVT, NPT, and production MD.")
# parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
# parser.add_argument('--gmxdir', '-g', dest='gmxdir', default='.', help='Directory containing GROMACS input files (ions.gro, topol.top, etc.)')
# parser.add_argument('--repn', '-n', dest='repn', type=int, default=None, help='Replicate number (used to create output directory)')
# parser.add_argument('--dna', dest='with_dna', action='store_true', help='Use if the system contains DNA. Adds DNA.bsc1 force field.')
# parser.add_argument('--phos', dest='with_phosphate', action='store_true', help='Use if the system contains phosphorylated residues.')
# parser.add_argument('--ntmpi', dest='ntmpi', type=int, default=0, help='Number of MPI ranks (default: auto-detected)')
# parser.add_argument('--ntomp', dest='ntomp', type=int, default=0, help='Number of OpenMP threads per rank (default: auto-detected)')
# parser.add_argument('--nb', dest='nb', type=str, default='gpu', help='Calculate non-bonded interactions on: auto, cpu, gpu')
# parser.add_argument('--pme', dest='pme', type=str, default='gpu', help='Perform PME calculations on: auto, cpu, gpu')
# parser.add_argument('--npme', dest='npme', type=int, default=-1, help='Number of separate ranks to be used for PME, -1 is guess (let GROMACS decide)')
# parser.add_argument('--ngpu', dest='ngpu', type=int, default=None, help='Number of GPUs to use (default: auto-detected)')
# parser.add_argument('--gputasks', dest='gputasks', type=str, default=None, help='Sequence of GPU IDs to use i.e. 0011 (default: auto-detected)')
# parser.add_argument('--em', dest='em', type=str, default=None, help='Energy minimization parameter file (default: minim.mdp)')
# parser.add_argument('--nvt', dest='nvt', type=str, default=None, help='NVT equilibration parameter file (default: nvt.mdp)')    
# parser.add_argument('--npt', dest='npt', type=str, default=None, help='NPT equilibration parameter file (default: npt.mdp)')
# parser.add_argument('--md', dest='md', type=str, default=None, help='Production MD parameter file (default: md.mdp)')

def main(logfile: Path, gmxdir: Path, outdir: Path, repn: int, 
         dna: bool, phos: bool, 
         ntmpi: int, ntomp: int, nb: str, pme: str, npme: int, 
         ngpu: int, gputasks: str, em: Path, nvt: Path, npt: Path, md: Path):

    ##### INITIALIZE LOG
    logger = setup_logger("md", logfile)

    ##### PARAMETER CHECKS #####

    logger.debug(f"Received positional argument gmxdir: {gmxdir} and replicate number: {repn}")

    # check if outdir exists, if not create it
    if not outdir.is_dir():
        outdir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Created output directory: {outdir}")

    # Check for checkpoint files
    checkpoint_files = list(outdir.glob(f"check_02_md_*.done"))
    if len(checkpoint_files) > 0:
        logger.warning("Checkpoint file found. Skipping MD steps: minimization, NVT, NPT and production MD.")
        logger.debug(f"To rerun these steps, please remove the checkpoint file(s): {outdir}/check_02_md_*.done")
        sys.exit(0)


    ##### RESOURCE SETUP #####
    logger.debug(f"Running {ntmpi} MPI ranks x {ntomp} threads per rank.")

    # Initialize empty flags
    flag_nb_gpu = []
    flag_pme_gpu = []
    flag_gpu_tasks = []
    GPUTASKS = ""

    # Assign flags based on resources
    if ngpu is not None:
        flag_nb_gpu = ["-nb", nb]
        flag_pme_gpu = ["-pme", pme, "-npme", str(npme)]
        if gputasks is not None:
            GPUTASKS = gputasks
            logger.debug(f"Using specified GPU list: {GPUTASKS}")
        else:
            logger.debug("No GPU task list specified. Generating GPU task list automatically.")
            for i in range(1, ntmpi+1):
                gpu = i % ngpu if ngpu > 0 else 0
                GPUTASKS += str(gpu)
        flag_gpu_tasks = ["-gputasks", GPUTASKS]
        logger.debug(f"Using {ngpu} GPU(s) with gputasks: {GPUTASKS}")
    else:
        logger.debug("No GPUs specified. Detecting available GPUs and assigning the default behavior.")


    ##### PARAMETER FILES SETUP #####
    # copy parameter files to outdir
    # parameter files
    param_mdp(logger, 'minim.mdp', outdir, em)
    param_mdp(logger, 'nvt.mdp', outdir, nvt)
    param_mdp(logger, 'npt.mdp', outdir, npt)
    param_mdp(logger, 'md.mdp', outdir, md)

    # check outdir contains minim.mdp, nvt.mdp, npt.mdp, md.mdp
    for mdp_file in ['minim.mdp', 'nvt.mdp', 'npt.mdp', 'md.mdp']:
        if not (outdir / mdp_file).is_file():
            logger.error(f"Parameter file {mdp_file} not found in output directory {outdir}.")
            sys.exit(1)

    logger.debug(f"Working directory: {Path.cwd()}, gmxdir: {gmxdir}, outdir: {outdir}")

    ##### MD STEPS #####

    ### 1. Energy Minimization
    if check_finished("em", outdir):
        logger.warning("em.log indicates minimization already completed and em.gro exists. Skipping energy minimization.")
    else:
        logger.debug("Starting energy minimization...")
        cmd = ["gmx", "grompp",
            "-f", outdir / "minim.mdp",
            "-c", gmxdir / "ions.gro",
            "-r", gmxdir / "ions.gro",
            "-p", gmxdir / "topol.top",
            "-o", outdir / "em.tpr",
            "-po", outdir / "em_out.mdp",
            "-quiet", "-maxwarn", "1"]
        subprocess.run(cmd, check=True, text=True)
        # log_process_output(process, "gmx grompp EM")

        # Calculate resource use
        em_ntmpi, em_ntomp = get_scaling(logger, outdir / "em.tpr", 'em')
        logger.debug(f"Running EM with: gmx mdrun -ntmpi {em_ntmpi} -ntomp {em_ntomp}")

        cmd = ["gmx", "mdrun", "-v", "-deffnm", outdir / "em", "-ntmpi", str(em_ntmpi), "-ntomp", str(em_ntomp), "-quiet"] + flag_nb_gpu
        subprocess.run(cmd, check=True, text=True)
        # log_process_output(process, "gmx mdrun EM")
        check_step(logger, "em", "Energy Minimization", outdir)

    # Check EM convergence and energy
    converged = False
    pot_energy = None
    fmax = None
    
    with open(outdir / "em.log", 'r') as f:
        content = f.read()
        # Extract final Potential Energy and Fmax
        pot_match = re.search(r'Potential Energy\s*=\s*([+-]?\d+\.?\d*[eE][+-]?\d+)', content)
        fmax_match = re.search(r'Maximum force\s*=\s*([+-]?\d+\.?\d*[eE][+-]?\d+)', content)
        conv_match = re.search(r'converged to Fmax', content)
        
        if pot_match:
            pot_energy = float(pot_match.group(1))
            pot_ok = pot_energy is not None and pot_energy < 0 and pot_energy < -1e5
        if fmax_match:
            fmax = float(fmax_match.group(1))
            fmax_ok = fmax is not None and fmax < 1e5  # Sanity check
        if conv_match:
            converged = True
        
    em_check = {
        'converged': converged,
        'pot_energy_ok': pot_ok,
        'pot_energy': pot_energy,
        'fmax': fmax,
        'success': converged and pot_ok and fmax_ok
    }
    if em_check['success']:
        logger.info("EM sanity check completed successfully.")
    else:
        logger.warning(f"EM sanity check failed. Consider tuning the EM parameters in a custom .mdp file.")
    logger.debug(f"Converged: {em_check['converged']}")
    logger.debug(f"PotE OK: {em_check['pot_energy_ok']} ({em_check['pot_energy']:.2e})")
    logger.debug(f"Fmax: {em_check['fmax']:.2e}")

    # Plot EM potential
    cmd = ["gmx", "energy", "-f", outdir / "em.edr", "-o", outdir / "em_potential.xvg", "-quiet"]
    subprocess.run(cmd, check=True, text=True, input="Potential\n0\n")

    x, y = np.loadtxt(outdir / "em_potential.xvg",comments=["@", "#"],unpack=True)
    plt.plot(x, y)
    plt.title("EM - Potential energy")
    plt.xlabel("Time (ps)")
    plt.ylabel("Potential energy (kJ/mol)")
    plt.savefig(outdir / "em_potential.png")
    plt.close()

    ### 2. Create groups
    group_index(gmxdir=outdir, logfile=logfile, logname='md', structurefile=outdir / "em.gro", dna=dna, phos=phos)

    ### 3. NVT Equilibration
    if check_finished("nvt", outdir):
        logger.warning("nvt.log indicates NVT equilibration already completed and nvt.gro exists. Skipping NVT equilibration.")
    else:
        logger.debug("Starting NVT equilibration...")
        cmd = ["gmx", "grompp",
                "-f", outdir / "nvt.mdp",
                "-c", outdir / "em.gro",
                "-r", outdir / "em.gro",
                "-p", gmxdir / "topol.top",
                "-o", outdir / "nvt.tpr",
                "-po", outdir / "nvt_out.mdp",
                "-n", outdir / "index.ndx", 
                "-quiet", "-maxwarn", "1"]
        subprocess.run(cmd, check=True, text=True)
        # log_process_output(process, "gmx grompp NVT")
        nvt_ntmpi, nvt_ntomp = get_scaling(logger, outdir / "nvt.tpr", 'nvt')
        subprocess.run(["gmx", "mdrun", "-v", "-deffnm", outdir / "nvt", "-ntmpi", str(nvt_ntmpi), "-ntomp", str(nvt_ntomp), "-quiet"] 
                    + flag_nb_gpu, check=True, text=True)
        # log_process_output(process, "gmx mdrun NVT")
        check_step(logger, "nvt", "NVT Equilibration", outdir)

    # Plot NVT Temperature
    cmd = ["gmx", "energy", "-f", outdir / "nvt.edr", "-o", outdir / "nvt_temp.xvg", "-quiet"]
    subprocess.run(cmd, check=True, text=True, input="Temperature\n0\n")

    x, y = np.loadtxt(outdir / "nvt_temp.xvg",comments=["@", "#"],unpack=True)
    final_seg_start = math.floor(max(x) / 3) * 2
    mean_final_seg = y[x >= final_seg_start].mean()
    std_final_seg = y[x >= final_seg_start].std()
    plt.plot(x, y)
    plt.axhline(y.mean(), color='b', ls='--', label = "Mean - full NVT")
    plt.axhline(mean_final_seg, color='r', ls='--', label = "Mean - last 1/3 time segment")
    plt.legend()
    plt.title("NVT Equilibration - Temperature")
    plt.xlabel("Time (ps)")
    plt.ylabel("Temperature (K)")
    plt.yticks(np.arange(round(min(y), -2), max(y), step=10))
    plt.savefig(outdir / "nvt_temp.png")
    plt.close()

    # Sanity check
    range_temp = range(290, 310)
    percent_out = np.sum((y[x >= final_seg_start] < range_temp.start) | (y[x >= final_seg_start] > range_temp.stop)) / len(y[x >= final_seg_start]) * 100
    logger.info(f"Temperature in the last 1/3 timepoints of NVT: average {mean_final_seg:.2f}, std {std_final_seg:.2f}.")
    logger.info(f"{percent_out:.2f}% of temperature values in the last 1/3 timepoints of NVT are outside of the default range ({range_temp.start}, {range_temp.stop}) (K).")
    if round(mean_final_seg) not in range_temp:
        logger.warning(f"Average temperature is outside of the default range ({range_temp.start}, {range_temp.stop}) (K), recommend for the default target temperature 300 K. " \
        "If this is your target temperature, consider tuning the NVT parameters in a custom .mdp file." \
        "If you are targeting a different temperature, please verify that the observed average temperature is close to your target.")
    else:
        logger.info(f"Average temperature is within the default range ({range_temp.start}, {range_temp.stop}) (K).")

    ### 4. NPT Equilibration
    if check_finished("npt", outdir):
        logger.warning("npt.log indicates NPT equilibration already completed and npt.gro exists. Skipping NPT equilibration.")
    else:
        logger.debug("Starting NPT equilibration...")
        cmd = ["gmx", "grompp",
            "-f", outdir / "npt.mdp",
            "-c", outdir / "nvt.gro",
            "-r", outdir / "nvt.gro",
            "-t", outdir / "nvt.cpt",
            "-p", gmxdir / "topol.top",
            "-o", outdir / "npt.tpr",
            "-po", outdir / "npt_out.mdp",
            "-n", outdir / "index.ndx",
            "-quiet", "-maxwarn", "1"]
        subprocess.run(cmd, check=True, text=True)
        # log_process_output(process, "gmx grompp NPT")
        npt_ntmpi, npt_ntomp = get_scaling(logger, outdir / "npt.tpr", 'npt')
        subprocess.run(["gmx", "mdrun", "-v", "-deffnm", outdir / "npt", "-ntmpi", str(npt_ntmpi), "-ntomp", str(npt_ntomp), "-quiet"] 
                    + flag_nb_gpu, check=True, text=True)
        # log_process_output(process, "gmx mdrun NPT")
        check_step(logger, "npt", "NPT Equilibration", outdir)

    # Plot NPT Pressure
    cmd = ["gmx", "energy", "-f", outdir / "npt.edr", "-o", outdir / "npt_pressure.xvg", "-quiet"]
    subprocess.run(cmd, check=True, text=True, input="Pressure\n0\n")

    x, y = np.loadtxt(outdir / "npt_pressure.xvg",comments=["@", "#"],unpack=True)
    final_seg_start = math.floor(max(x) / 3) * 2
    mean_final_seg = y[x >= final_seg_start].mean()
    std_final_seg = y[x >= final_seg_start].std()
    plt.plot(x, y)
    plt.axhline(y.mean(), color='b', ls='--', label = "Mean - full NPT")
    plt.axhline(mean_final_seg,color='r', ls='--', label = "Mean - last 1/3 time segment")
    plt.legend()
    plt.title("NPT Equilibration - Pressure")
    plt.xlabel("Time (ps)")
    plt.ylabel("Pressure (bar)")
    plt.yticks(np.arange(round(min(y), -3), max(y) , step=200))
    plt.savefig(outdir / "npt_pressure.png")
    plt.close()

    # Sanity check
    range_pres = range(-11, 11)
    percent_out = np.sum((y[x >= final_seg_start] < range_pres.start) | (y[x >= final_seg_start] > range_pres.stop)) / len(y[x >= final_seg_start]) * 100
    logger.info(f"Pressure in the last 1/3 timepoints of NPT: average {mean_final_seg:.3f}, std {std_final_seg:.3f}.")
    logger.info(f"{percent_out:.2f}% of pressure values in the last 1/3 timepoints of NPT are outside of the default range ({range_pres.start}, {range_pres.stop}) (bar).")
    if round(mean_final_seg) not in range_pres:
        logger.warning(f"Average pressure is outside of the default range ({range_pres.start}, {range_pres.stop}) (bar). " \
        "If this is your target pressure, consider tuning the NPT parameters in a custom .mdp file. " \
        "If you are targeting a different pressure, please verify that the observed average pressure is close to your target.")
    else:
        logger.info(f"Average pressure is within the default range ({range_pres.start}, {range_pres.stop}) (bar).")

    # Plot NPT Density
    cmd = ["gmx", "energy", "-f", outdir / "npt.edr", "-o", outdir / "npt_density.xvg", "-quiet"]
    subprocess.run(cmd, check=True, text=True, input="Density\n0\n")

    x, y = np.loadtxt(outdir / "npt_density.xvg",comments=["@", "#"],unpack=True)
    final_seg_start = math.floor(max(x) / 3) * 2
    mean_final_seg = y[x >= final_seg_start].mean()
    std_final_seg = y[x >= final_seg_start].std()
    plt.plot(x, y)
    plt.axhline(y.mean(), color='b', ls='--', label = "Mean - full NPT")
    plt.axhline(mean_final_seg,color='r', ls='--', label = "Mean - last 1/3 time segment")
    plt.legend()
    plt.title("NPT Equilibration - Density")
    plt.xlabel("Time (ps)")
    plt.ylabel("Density (kg/m\u00b3)")
    plt.savefig(outdir / "npt_density.png")
    plt.close()

    # Sanity check
    logger.info(f"Density in the last 1/3 timepoints of NPT: average {mean_final_seg:.3f}, std {std_final_seg:.3f}.")
    if round(std_final_seg) not in range(-5, 5):
        logger.warning("Standard deviation of density in the last 1/3 timepoints of NPT is outside of the default range (-5, 5) (kg/m^3). " \
        "This might indicate that the system did not stabilize. "  \
        "Check the density plot and consider tuning the NPT parameters in a custom .mdp file.")
    else:
        logger.info("Standard deviation of density is within the default range (-5, 5) (kg/m^3).")


    ### 5. Production MD
    if check_finished("md", outdir):
        logger.warning("md.log indicates production MD already completed and md.gro exists. Skipping production MD.")
    else:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        logfile = os.path.join(outdir, f"md_{timestamp}.log")
        if os.path.isfile(f"{outdir}/md.cpt"):
            logger.debug(f"Resuming production MD from checkpoint {outdir}/md.cpt...")
            cmd = ["gmx", "mdrun", 
                "-ntmpi", str(ntmpi), 
                "-ntomp", str(ntomp)] + flag_nb_gpu + flag_pme_gpu + flag_gpu_tasks + [
                    "-cpi", outdir / "md.cpt", "-append", "-v", "-deffnm", outdir / "md", "-quiet"]
            subprocess.run(cmd, check=True, text=True)
            # log_process_output(process, "gmx mdrun MD (from checkpoint)")
        else:
            logger.debug("Starting production MD...")
            cmd = ["gmx", "grompp",
                "-f", outdir / "md.mdp",
                "-c", outdir / "npt.gro",
                "-r", outdir / "npt.gro",
                "-t", outdir / "npt.cpt",
                "-p", gmxdir / "topol.top",
                "-o", outdir / "md.tpr",
                    "-po", outdir / "md_out.mdp",
                "-n", outdir / "index.ndx",
                "-quiet", "-v", "-maxwarn", "1"]
            subprocess.run(cmd, check=True, text=True)
            # log_process_output(process, "gmx grompp MD")
            if ntmpi == 0 or ntomp == 0:
                logger.warning("ntmpi or ntomp is set to 0. Calculating safe defaults based on system size and available resources.")
                ntmpi, ntomp = get_scaling(logger, outdir / "md.tpr", 'md')
                logger.debug(f"Calculated scaling for production MD: ntmpi={ntmpi}, ntomp={ntomp}")
            cmd = ["gmx", "mdrun", "-ntmpi", str(ntmpi), "-ntomp", str(ntomp)
                ] + flag_nb_gpu + flag_pme_gpu + flag_gpu_tasks + ["-v", "-deffnm", outdir / "md", "-quiet"]
            logger.debug(f"Production MD command: {cmd}")
            subprocess.run(cmd, check=True, text=True)
            # log_process_output(process, "gmx mdrun MD")
        check_step(logger, "md", "Production MD", outdir)

    # Create checkpoint file
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    checkpoint = outdir / f"check_02_md_{timestamp}.done"
    checkpoint.parent.mkdir(exist_ok=True, parents=True)
    with checkpoint.open('w') as f:
        f.write("MD completed successfully.\n")
        logger.info("MD steps completed successfully. Checkpoint created: " + str(checkpoint))
        return


if __name__ == "__main__":
    main()