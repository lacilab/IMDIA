#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import multiprocessing
import subprocess
import shutil
from sympy import isprime
from pathlib import Path

##### HELPER FUNCTIONS FOR MD #####

def get_total_cores():
    if 'SLURM_CPUS_PER_TASK' in os.environ:
        return int(os.environ['SLURM_CPUS_PER_TASK'])
    elif 'SLURM_NTASKS' in os.environ:
        cpus_per_task = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
        return int(os.environ['SLURM_NTASKS']) * cpus_per_task
    
    # Fallback: physical cores
    return multiprocessing.cpu_count() // 2

def param_mdp(logger, filename: str, outdir: Path, custom_path: Path = None):
    """
    Set parameter files. args: <custom_path_var> <filename>
    If -em, --nvt, --npt, --md options are provided, copy those files instead of the default ones to $outdir. 
    Else, copy the default files
    """
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

def LastNlines(fname, N):
    """
    Helper to read last N lines of the file 
    (adapted from GeeksforGeeks)
    """
    lines = []
    with open(fname) as file:
        for line in (file.readlines()[-N:]):
            lines.append(line)
    return lines

def check_finished(step_short: str, outdir: Path):
    """
    Check .gro and .log files to check if a mdrun completed successfully.

    Arguments:
        stepshort : Short name of the step, e.g. em, nvt, npt, md.
        outdir    : Path to the directory with the MD results (e.g. rep1/gromacs).
    """
    if (outdir / f"{step_short}.gro").is_file() and (outdir / f"{step_short}.log").is_file() and "Finished mdrun" in " ".join(LastNlines(outdir / f"{step_short}.log", 20)):
        return True
    else:
        return False
    
def log_process_output(process, step_name, logger):
    if process.returncode != 0:
        logger.error(f"{step_name} error:\n {process.stderr}")
    elif process.stdout:
        logger.debug(f"{step_name} output:\n {process.stdout}")   


def check_step(logger, step_short: str, step_name: str, outdir: Path):
    """
    Log step completion. 

    Arguments:
        logger        : Logger object.
        stepshort : Short name of the step, e.g. em, nvt, npt, md.
        step_name : Long name of the step, e.g. energy minimization.
        outdir    : Path to the directory with the MD results (e.g. rep1/gromacs).
    """
    if check_finished(step_short, outdir):
        logger.info(f"{step_name} completed.")
    else:
        logger.error(f"{step_name} failed. Check {outdir}/{step_short}.log for details.")
        sys.exit(1)

def get_grid_size(logger, stage:str, dump): 
    """ 
    Get fourier grid dimensions and use them to compute the box size.

    Arguments:
        logger  : Logger object.
        dump    : gmx dump output from tpr file.
    """   
    # Parse separate fourier variables
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
    """
    Extract the number of atoms from a GROMACS dump output.

    Arguments:
        logger  : Logger object.
        dump    : gmx dump output from tpr file.

    Returns:
        natoms  : Number of atoms parsed from the dump. Defaults to 10000 if not found.
    """
    natoms_match = re.search(r'natoms\s*=\s*(\d+)', dump)
    natoms = int(natoms_match.group(1)) if natoms_match else 10000
    logger.debug(f"Parsed number of atoms: {natoms}")
    return natoms

def get_mincellsize(logger, logfile: Path):
    """
    Extract minimum cell size from EM log.
    
    Arguments:
        logger  : Logger object.
        logfile : mdrun log file from an MD step (usually EM)
    """

    if not logfile.is_file():
        logger.warning(f"Log file '{logfile}' not found. Cannot extract minimum cell size for scaling decisions.")
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
    """
    Determine optimal parallelization settings (ntmpi and ntomp) for a GROMACS simulation stage
    based on available computational resources and system properties.

    This function parses the input .tpr file using `gmx dump` to extract information such as
    the simulation box size and number of atoms, and optionally uses the minimum cell size
    from the energy minimization log. It then computes a safe and efficient combination of
    MPI ranks (ntmpi) and OpenMP threads (ntomp) to balance performance and avoid domain
    decomposition errors.

    Arguments:
        logger       : Logger object for debug output.
        tpr_file     : Path to the GROMACS .tpr file.
        stage        : Simulation stage (e.g., "em", "nvt", "npt", "md").
        total_cores  : Total number of CPU cores available (default: auto-detected).

    Returns:
        ntmpi        : Number of MPI ranks.
        ntomp        : Number of OpenMP threads per rank.
    """
    
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