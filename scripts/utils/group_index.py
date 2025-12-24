#!/usr/bin/env python

import os
import re
import subprocess
from scripts.utils.helpers import setup_logger


def group_index(gmxdir: str = None, logfile: str = None, logname: str = 'log', structurefile: str = None, dna: bool = False, phos: bool = False):
    
    def log_process_output(process, step_name):
        if process.returncode != 0:
            logger.error(f"{step_name} error:\n {process.stderr}")
        elif process.stdout:
            logger.debug(f"{step_name} output:\n {process.stdout}")   
    
    ##### INITIALIZE LOG
    logger = setup_logger(logname, logfile)

    # Construct group name
    group = "Protein"
    if dna:
        group += "_DNA"
    if phos:
        group += "_P"

    # Substitute group name in nvt.mdp, npt.mdp, md.mdp
    if group != 'Protein':
        for mdp_file in ['nvt.mdp', 'npt.mdp', 'md.mdp']:
            mdp_path = os.path.join(gmxdir, mdp_file)
            with open(mdp_path, 'r') as f:
                contents = f.read()
            contents = contents.replace("Protein", group)
            with open(mdp_path, 'w') as f:
                f.write(contents)

    # Check if an index file with the target groups already exists
    ndx_ready = False
    target_groups = [group, 'Chain_A', 'Chain_B', 'Water_Ions']
    if os.path.isfile(f"{gmxdir}/index.ndx"):
        logger.debug(f'Index file {gmxdir}/index.ndx already exists. Checking for target groups {target_groups}')
        with open(f"{gmxdir}/index.ndx") as f:
            ndx_groups = []
            for line in f:
                line = line.strip()
                if line.startswith("[") and line.endswith("]"):
                    ndx_groups.append(line[1:-1].strip())
            logger.debug(f'Groups found in existing index file: {ndx_groups}')
        if set(target_groups).issubset(ndx_groups):
            logger.warning(f'The index file {gmxdir}/index.ndx already exists and contains the target groups {target_groups}. Skipping index creation step.')
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
        group_water = groups.index("Water")
        group_na = groups.index("Na")
        group_cl = groups.index("Cl")
        # if dna:
        #     group_dna = groups.index("DNA")
        if phos:
            group_phosphate = groups.index("PTR")

        # Remove index file
        os.remove(f"{gmxdir}/temp.ndx")
        # Run again to create index file with new groups
        # groups_input = f"""
        # chain A
        # name {group_n} Chain_A
        # chain B
        # name {group_n + 1} Chain_B
        # """
        groups_input = ""
        if dna and phos:
            groups_input += f"""
            r DA | r DC | r DT | r DG | r DC3 | r DC5 | r DG3 | r DG5
            name {group_n} DNA
            {group_protein} | {group_n} | {group_phosphate}
            name {group_n + 1} Protein_DNA_P
            {group_protein} | {group_phosphate}
            name {group_n + 2} Protein_P
            {group_water} | {group_na} | {group_cl}
            name {group_n + 3} Water_Ions
            quit
            """
        elif dna and not phos:    
            groups_input += f"""
            r DA | r DC | r DT | r DG | r DC3 | r DC5 | r DG3 | r DG5
            name {group_n} DNA
            {group_protein} | {group_n}
            name {group_n + 1} Protein_DNA
            {group_water} | {group_na} | {group_cl}
            name {group_n + 2} Water_Ions
            quit
            """
        else:
            groups_input += f"""
            {group_water} | {group_na} | {group_cl}
            name {group_n} Water_Ions
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
        ndx_ready = True

    return group
    
