from pathlib import Path
import subprocess
import sys
import os
import re
import datetime
import math
import numpy as np
import matplotlib.pyplot as plt
from scripts.utils.md_helpers import param_mdp, check_finished, check_step, get_scaling
from scripts.utils.helpers import setup_logger
from scripts.utils.group_index import group_index

def setup(logfile: Path, gmxdir: Path, outdir: Path, repn: int,
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
    param_mdp(logger, 'em.mdp', outdir, em)
    param_mdp(logger, 'nvt.mdp', outdir, nvt)
    param_mdp(logger, 'npt.mdp', outdir, npt)
    param_mdp(logger, 'md.mdp', outdir, md)

    # check outdir contains em.mdp, nvt.mdp, npt.mdp, md.mdp
    for mdp_file in ['em.mdp', 'nvt.mdp', 'npt.mdp', 'md.mdp']:
        if not (outdir / mdp_file).is_file():
            logger.error(f"Parameter file {mdp_file} not found in output directory {outdir}.")
            sys.exit(1)

    logger.debug(f"Working directory: {Path.cwd()}, gmxdir: {gmxdir}, outdir: {outdir}")
    return  flag_nb_gpu, flag_pme_gpu, flag_gpu_tasks

##### MD STEPS #####
### 1. Energy Minimization
def minimization(logfile: Path, gmxdir: Path, outdir: Path, flag_nb_gpu: list):
    
    ##### INITIALIZE LOG
    logger = setup_logger("md", logfile)
    
    if check_finished("em", outdir):
        logger.warning("em.log indicates minimization already completed and em.gro exists. Skipping energy minimization.")
    else:
        logger.debug("Starting energy minimization...")
        cmd = ["gmx", "grompp",
            "-f", outdir / "em.mdp",
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
    if (outdir / "em.edr").is_file() and not (outdir / "em_potential.png").is_file():
        cmd = ["gmx", "energy", "-f", outdir / "em.edr", "-o", outdir / "em_potential.xvg", "-quiet"]
        subprocess.run(cmd, check=True, text=True, input="Potential\n0\n")
    if (outdir / "em_potential.png").is_file() and not (outdir / "em_potential.png").is_file():
        x, y = np.loadtxt(outdir / "em_potential.xvg",comments=["@", "#"],unpack=True)
        plt.plot(x, y)
        plt.title("EM - Potential energy")
        plt.xlabel("Time (ps)")
        plt.ylabel("Potential energy (kJ/mol)")
        plt.savefig(outdir / "em_potential.png")
        plt.close()



def equilibration(logfile: Path, gmxdir: Path, outdir: Path, flag_nb_gpu: list):
    
    ##### INITIALIZE LOG
    logger = setup_logger("md", logfile)

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
    if (outdir / "nvt.edr").is_file() and not (outdir / "nvt_temp.xvg").is_file():
        cmd = ["gmx", "energy", "-f", outdir / "nvt.edr", "-o", outdir / "nvt_temp.xvg", "-quiet"]
        subprocess.run(cmd, check=True, text=True, input="Temperature\n0\n")

    x, y = np.loadtxt(outdir / "nvt_temp.xvg",comments=["@", "#"],unpack=True)
    final_seg_start = math.floor(max(x) / 3) * 2
    mean_final_seg = y[x >= final_seg_start].mean()
    std_final_seg = y[x >= final_seg_start].std()

    if not (outdir / "nvt_temp.png").is_file():
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
    if (outdir / "npt.edr").is_file() and not (outdir / "npt_pressure.xvg").is_file():
        cmd = ["gmx", "energy", "-f", outdir / "npt.edr", "-o", outdir / "npt_pressure.xvg", "-quiet"]
        subprocess.run(cmd, check=True, text=True, input="Pressure\n0\n")

    x, y = np.loadtxt(outdir / "npt_pressure.xvg",comments=["@", "#"],unpack=True)
    final_seg_start = math.floor(max(x) / 3) * 2
    mean_final_seg = y[x >= final_seg_start].mean()
    std_final_seg = y[x >= final_seg_start].std()

    if not (outdir / "npt_pressure.png").is_file():
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
    if not (outdir / "npt_density.xvg").is_file():
        cmd = ["gmx", "energy", "-f", outdir / "npt.edr", "-o", outdir / "npt_density.xvg", "-quiet"]
        subprocess.run(cmd, check=True, text=True, input="Density\n0\n")

    x, y = np.loadtxt(outdir / "npt_density.xvg",comments=["@", "#"],unpack=True)
    final_seg_start = math.floor(max(x) / 3) * 2
    mean_final_seg = y[x >= final_seg_start].mean()
    std_final_seg = y[x >= final_seg_start].std()

    if not (outdir / "npt_density.png").is_file():
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



def production(logfile: Path, gmxdir: Path, outdir: Path, ntmpi: int, ntomp: int, flag_nb_gpu: list, flag_pme_gpu: list, flag_gpu_tasks: list):
    
    ##### INITIALIZE LOG
    logger = setup_logger("md", logfile)

    ### 5. Production MD
    if check_finished("md", outdir):
        logger.warning("md.log indicates production MD already completed and md.gro exists. Skipping production MD.")
    else:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        logfile = os.path.join(outdir, f"md_{timestamp}.log")
        if os.path.isfile(f"{outdir}/md.cpt"):
            if ntmpi == 0 or ntomp == 0:
                logger.warning("ntmpi or ntomp is set to 0. Calculating safe defaults based on system size and available resources.")
                ntmpi, ntomp = get_scaling(logger, outdir / "md.tpr", 'md')
                logger.debug(f"Calculated scaling for production MD: ntmpi={ntmpi}, ntomp={ntomp}")
            logger.debug(f"Resuming production MD from checkpoint {outdir}/md.cpt...")
            cmd = ["gmx", "mdrun", 
                "-ntmpi", str(ntmpi), 
                "-ntomp", str(ntomp)] + flag_nb_gpu + flag_pme_gpu + flag_gpu_tasks + [
                    "-cpi", outdir / "md.cpt", "-append", "-v", "-deffnm", outdir / "md", "-quiet"]
            subprocess.run(cmd, check=True, text=True)
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
