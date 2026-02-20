#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import subprocess
from scripts.utils.helpers import setup_logger

def main(logfile: Path, datadir: Path, configfile: Path, cores: int = None, rules: list[str] = None, unlock: bool = False):
    # parser = argparse.ArgumentParser(description="Analyse MD simulation results.")
    # parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
    # parser.add_argument('--datadir', '-d', dest='datadir', default=None, help='Path to the data directory.')
    # parser.add_argument('--configfile', '-c', dest='configfile', type=str, default=None, 
    #                     help='Path to the YAML configuration file to be used in snakefile.')
    # parser.add_argument('--cores', '-n', dest='cores', default="1", 
    #                     help='Number of cores to use for analysis. Default: all available CPUs.')
    # parser.add_argument('--rules', dest='rules', nargs='+', default=None, 
    #                     help='Snakefile rules to be executed.')
    # parser.add_argument('--ffs', dest='ffs', type=str, default=None,
    #                     help='Force fields used in GBSA.')
    # parser.add_argument('--unlock', dest='unlock', action='store_true', help='Unlock snakemake workflow if it is locked.')

    # Initialize logger
    logger = setup_logger("analysis", logfile)
    logger.debug(f"Starting analysis with snakemake using config file: {configfile}.")
    
    if rules is not None:
        logger.debug(f"Executing specific rules: {rules}")
    else:
        logger.debug("No rules were specified. Executing all rules in the Snakefile.")
        
    config_args = [f"logfile={str(logfile)}", f"datadir={str(datadir)}"]

    if rules:
        config_args.append(f"rule_groups={' '.join(rules)}")

    if unlock:
        logger.debug("Unlocking snakemake workflow.")
        cmd = [
            "snakemake",
            "--snakefile", "scripts/analysis04/Snakefile",
            "--cores", "1",
            "--configfile", str(configfile),
            "--config",
            *config_args,
            "--unlock"
        ]
    else:
        cmd = [
            "snakemake",
            "--snakefile", "scripts/analysis04/Snakefile",
            "--cores", str(cores) if cores is not None else "all",
            "--configfile", str(configfile),
            "--config",
            *config_args
        ]

    # process = subprocess.run(cmd, text=True, check=True, capture_output=True)
    # if process.returncode != 0:
    #     logger.error(f"snakemake error:\n {process.stderr}")
    # elif process.stdout:
    #     logger.debug(f"snakemake output:\n {process.stdout}")  

    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()