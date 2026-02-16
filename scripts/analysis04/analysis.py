#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import subprocess
import sys
from scripts.utils.helpers import setup_logger

def main():
    parser = argparse.ArgumentParser(description="Analyse MD simulation results.")
    parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
    parser.add_argument('--datadir', '-d', dest='datadir', default=None, help='Path to the data directory.')
    parser.add_argument('--configfile', '-c', dest='configfile', type=str, default=None, 
                        help='Path to the YAML configuration file to be used in snakefile.')
    parser.add_argument('--cores', '-n', dest='cores', default="1", 
                        help='Number of cores to use for analysis. Default: all available CPUs.')
    parser.add_argument('--rules', dest='rules', nargs='+', default=None, 
                        help='Snakefile rules to be executed.')
    parser.add_argument('--ffs', dest='ffs', type=str, default=None,
                        help='Force fields used in GBSA.')
    parser.add_argument('--unlock', dest='unlock', action='store_true', help='Unlock snakemake workflow if it is locked.')

    # print help if required arguments are not provided
    if len(sys.argv) < 5:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Initialize logger
    logger = setup_logger("analysis", args.logfile)
    logger.debug(f"Starting analysis with snakemake using config file: {args.configfile}.")
    
    if args.rules is not None:
        logger.debug(f"Executing specific rules: {args.rules}")
    else:
        logger.debug("No rules were specified. Executing all rules in the Snakefile.")
        
    config_args = [f"logfile={args.logfile}", f"datadir={args.datadir}"]

    if args.rules is not None:
        config_args.append(f"rule_groups={' '.join(args.rules)}")

    if args.ffs is not None:
        config_args.append(f"ffs={' '.join(args.ffs)}")

    if args.unlock:
        logger.debug("Unlocking snakemake workflow.")
        cmd = [
            "snakemake",
            "--snakefile", "scripts/analysis04/Snakefile",
            "--cores", "1",
            "--configfile", args.configfile,
            "--config",
            *config_args,
            "--unlock"
        ]
    else:
        cmd = [
            "snakemake",
            "--snakefile", "scripts/analysis04/Snakefile",
            "--cores", str(args.cores) if args.cores else "all",
            "--configfile", args.configfile,
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