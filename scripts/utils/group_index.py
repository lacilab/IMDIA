#!/usr/bin/env python

import os
import re
import subprocess
import sys
from pathlib import Path
from scripts.utils.helpers import setup_logger


def group_index(gmxdir: Path = None, logfile: Path = None, logname: str = 'log', structurefile: Path = None, dna: bool = False, phos: bool = False):
    
    def log_process_output(process, step_name):
        if process.returncode != 0:
            logger.error(f"{step_name} error:\n {process.stderr}")
        elif process.stdout:
            logger.debug(f"{step_name} output:\n {process.stdout}")   
    
    ##### INITIALIZE LOG
    logger = setup_logger(logname, logfile)

    METALS = ['ZN','MG','FE2','FE3','CA','MN','CU','NI','CO','CD','PB','HG','LI']

    # Construct group name
    group = "Protein"
    if dna:
        group += "_DNA"
    if phos:
        group += "_P"

    # In the MD step, substitute group name in nvt.mdp, npt.mdp, md.mdp
    if logname == 'md' and group != 'Protein':
        for mdp_file in ['nvt.mdp', 'npt.mdp', 'md.mdp']:
            mdp_path = gmxdir / mdp_file
            with mdp_path.open('r') as f:
                contents = f.read()
            contents = contents.replace("Protein", group)
            with mdp_path.open('w') as f:
                f.write(contents)

    # Check if an index file with the target groups already exists
    target_groups = [group, 'Water_Ions']
    index_ori = gmxdir / "index.ndx"
    if index_ori.is_file():
        logger.debug(f'Index file {gmxdir}/index.ndx already exists. Checking for target groups {target_groups}')
        with index_ori.open('r') as f:
            ndx_groups = []
            for line in f:
                line = line.strip()
                if line.startswith("[") and line.endswith("]"):
                    ndx_groups.append(line[1:-1].strip())
            logger.debug(f'Groups found in existing index file: {ndx_groups}')
        if set(target_groups).issubset(ndx_groups):
            logger.warning(f'The index file {gmxdir}/index.ndx already exists and contains the target groups {target_groups}. Skipping index creation step.')
            return group

    index_temp = gmxdir / "temp.ndx"
    # Get groups in the index file
    cmd = [
        "gmx", "make_ndx",
        "-f", structurefile,
        "-o", index_temp,
        "-quiet"
    ]

    # Run gmx make_ndx and auto-quit with "q"
    result = subprocess.run(
        cmd,
        input="q\n",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    index_output = result.stdout

    # Regex to extract group names
    group_pattern = re.compile(r"^\s*\d+\s+([A-Za-z0-9_+\-]+)\s*:", re.MULTILINE)
    groups = group_pattern.findall(index_output)
    group_n = len(groups)
    logger.debug(f"Found {group_n} groups in index file {gmxdir}/temp.ndx: {groups}")

    group_protein = groups.index("Protein")
    group_water = groups.index("Water")
    group_na = groups.index("Na")
    group_cl = groups.index("Cl")

    if phos:
        group_phosphate = groups.index("PTR")

    # Remove index file
    os.remove(index_temp)

    groups_input = ""
    if dna and phos:
        groups_input += f"""
        r DA | r DC | r DT | r DG | r DC3 | r DC5 | r DG3 | r DG5
        name {group_n} DNA
        {"r " + " | r ".join(METALS)}
        name {group_n + 1} Metals
        {group_protein} | {group_n} | {group_phosphate} | {group_n + 1}
        name {group_n + 2} Protein_DNA_P
        {group_protein} | {group_phosphate}
        name {group_n + 3} Protein_P
        {group_water} | {group_na} | {group_cl}
        name {group_n + 4} Water_Ions
        quit
        """
    elif dna and not phos:    
        groups_input += f"""
        r DA | r DC | r DT | r DG | r DC3 | r DC5 | r DG3 | r DG5
        name {group_n} DNA
        {"r " + " | r ".join(METALS)}
        name {group_n + 1} Metals
        {group_protein} | {group_n} | {group_n + 1}
        name {group_n + 2} Protein_DNA
        {group_water} | {group_na} | {group_cl}
        name {group_n + 3} Water_Ions
        quit
        """
    else:
        groups_input += f"""
        {"r " + " | r ".join(METALS)}
        name {group_n} Metals
        {group_protein} | {group_n}
        name {group_n + 1} Protein
        {group_water} | {group_na} | {group_cl}
        name {group_n + 2} Water_Ions
        quit
        """

    cmd = [ "gmx", "make_ndx", 
        "-f", structurefile,
        "-o", f"{gmxdir}/temp_groups.ndx",
        "-quiet"
        ]
    process = subprocess.run(cmd, input=groups_input, text=True, capture_output=True, check=True)
    log_process_output(process, "gmx make_ndx")
    logger.info("Index file created with new groups")
    # Move temp_groups.ndx to final index.ndx
    os.replace(f"{gmxdir}/temp_groups.ndx", f"{gmxdir}/index.ndx")
    logger.info(f"Index file with target groups saved as {gmxdir}/index.ndx")
    return group


def group_index_interactions(gmxdir: str = None, logfile: str = None, logname: str = 'log', structurefile: str = None, indexfile: str = None,
                             group1: str = None, group2: str = None, dna: bool = False, phos: bool = False):
    
    def log_process_output(process, step_name):
        if process.returncode != 0:
            logger.error(f"{step_name} error:\n {process.stderr}")
        elif process.stdout:
            logger.debug(f"{step_name} output:\n {process.stdout}")   
    
    ##### INITIALIZE LOG
    logger = setup_logger(logname, logfile)

    # Check that group1 and group2 are provided
    if group1 is None or group2 is None:
        logger.error("Both group1 and group2 must be specified for interaction analysis.")
        sys.exit(1)

    # Check if an index file with the target groups already exists
    ndx_ready = False
    ndx_names = []
    target_groups = [group1, group2]
    for i, group in enumerate([group1, group2]):
            if re.match(r"^[A-Za-z]$", group):
                ndx_names.append(f"Chain_{group}")
            elif group == "DNA" and dna:
                ndx_names.append("DNA")
            elif group == "Protein" and phos:
                ndx_names.append("Protein_P")
            else: # Protein group already predifined   
                ndx_names.append("Protein")
    if os.path.isfile(indexfile):
        logger.debug(f'Index file {indexfile} already exists. Checking for target groups {ndx_names}')
        with open(indexfile) as f:
            ndx_groups = []
            for line in f:
                line = line.strip()
                if line.startswith("[") and line.endswith("]"):
                    ndx_groups.append(line[1:-1].strip())
            logger.debug(f'Groups found in existing index file: {ndx_groups}')
        if set(ndx_names).issubset(ndx_groups):
            logger.warning(f'The index file {indexfile} already exists and contains the target groups {ndx_names}. Skipping index creation step.')
            ndx_names = target_groups.copy()
            ndx_ready = True

    if not ndx_ready:
        # Get groups in the index file
        cmd = [
            "gmx", "make_ndx",
            "-f", structurefile,
            "-o", f"{gmxdir}/temp.ndx",
            "-quiet"
        ]

        # Run gmx make_ndx and auto-quit with "q"
        result = subprocess.run(
            cmd,
            input="q\n",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        index_output = result.stdout

        # Regex to extract group names
        group_pattern = re.compile(r"^\s*\d+\s+([A-Za-z0-9_+\-]+)\s*:", re.MULTILINE)
        groups = group_pattern.findall(index_output)
        group_n = len(groups)
        logger.debug(f"Found {group_n} groups in index file {gmxdir}/temp.ndx: {groups}")

        group_protein = groups.index("Protein")
        if phos:
            group_phosphate = groups.index("PTR")

        # Remove index file
        os.remove(f"{gmxdir}/temp.ndx")
        # Run again to create index file with new groups
        groups_input = ""
        # if group1 and group2 are single letters (representing chains), map them to the corresponding chains
        for i, group in enumerate([group1, group2]):
            if re.match(r"^[A-Za-z]$", group):
                logger.debug(f"Creating group for chain {group} as number {group_n + i}")
                groups_input += f"""chain {group}
                name {group_n + i} Chain_{group}
"""
            elif group == "DNA" and dna:
                groups_input += f"""
                            r DA | r DC | r DT | r DG | r DC3 | r DC5 | r DG3 | r DG5
                            name {group_n + i} DNA
                            """
            elif group == "Protein" and phos:
                groups_input += f"""
                            {group_protein} | {group_phosphate}
                            name {group_n + i} Protein_P
                            """
                 
        
        groups_input += """
            quit
            """

        cmd = [ "gmx", "make_ndx", 
            "-f", structurefile,
            "-o", f"{gmxdir}/temp_groups.ndx",
            "-quiet"
            ]
        subprocess.run(cmd, input=groups_input, text=True, check=True)
        # log_process_output(process, "gmx make_ndx")
        logger.info(f"Index file created with new groups")

        # Move temp_groups.ndx to final index.ndx
        os.replace(f"{gmxdir}/temp_groups.ndx", indexfile)
        logger.info(f"Index file with target groups saved as {indexfile}")
        ndx_ready = True

        # Check that ndx_names contains two names
        if len(ndx_names) != 2:
            logger.error(f"Error in creating list of index group names. Expected 2 names, got {len(ndx_names)}.")
            sys.exit(1)

    return ndx_names
