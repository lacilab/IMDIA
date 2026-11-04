#!/usr/bin/env python

import os
import re
import subprocess
import sys
from pathlib import Path
import shutil
import MDAnalysis as mda
from scripts.utils.helpers import setup_logger

def log_process_output(logger, process, step_name):
    if process.returncode != 0:
        logger.error(f"{step_name} error:\n {process.stderr}")
    elif process.stdout:
        logger.debug(f"{step_name} output:\n {process.stdout}")   

def get_name_main(logger, structurefile: Path = None, dna: bool = False, rna: bool = False, phos: bool = False):
    """
    Build the main group name based on the presence of DNA, RNA, phosphate groups, and metal ions in the structure file.
    
    Arguments:
        logger        : Logger object.
        structurefile : Path to PDB structure file.
        dna           : True if DNA is present in the system.
        rna           : True if RNA is present in the system.
        phos          : True if the system contains phosphorylated residues.
    """
    
    # Identify metals from list present in the structure file
    METALS = ['ZN','MG','FE2','FE3','CA','MN','CU','NI','CO','CD','PB','HG','LI']

    u = mda.Universe(structurefile)
    residue_names_in_structure = sorted(set(u.atoms.resnames))

    present_metals = [metal for metal in METALS if metal in residue_names_in_structure]
    if present_metals:
        logger.debug(f"Metal ions found in structure: {present_metals}")
    else:
        logger.debug("No metal ions found in structure.")

    # Identify DNA and RNA residues from the names present in the structure file
    # (AMBER names: DA5/DA/DA3/DAN, A5/A/A3/AN; GROMACS/CHARMM names: RA5/RA/RA3/RAN).
    DNA_RESIDUES = ['DA','DC','DT','DG','DC3','DC5','DG3','DG5', "DA3", "DA5", "DT3", "DT5", "DU", "DU3", "DU5",
                    "DAN", "DCN", "DGN", "DTN"]
    RNA_RESIDUES = [f"{prefix}{base}{suffix}" for prefix in ("", "R") for base in "ACGU"
                    for suffix in ("", "5", "3", "N")]
    present_dna = [nt for nt in DNA_RESIDUES if nt in residue_names_in_structure]
    present_rna = [nt for nt in RNA_RESIDUES if nt in residue_names_in_structure]
    if present_dna or present_rna:
        logger.debug(f"Nucleotides found in structure: DNA {present_dna}, RNA {present_rna}")
    else:
        logger.debug("No nucleotides found in structure.")

    # Construct main group name
    group_main = "Protein"
    group_protein = "Protein"
    if dna:
        group_main += "_DNA"
    if rna:
        group_main += "_RNA"
    if phos:
        group_main += "_P"
        group_protein += "_P"
    if len(present_metals) > 0:
        group_main += "_Metals"
        group_protein += "_Metals"

    return group_main, group_protein, present_metals, present_dna, present_rna

def targets_in_index(logger, gmxdir: Path, target_groups: list[str]):
    """
    Check if an index file with the target groups already exists.
    
    Arguments:
        logger        : Logger object.
        gmxdir        : Path to the directory with the MD results (e.g. rep1/gromacs).
        target_groups : List of target groups.
    """
    targets_comp = target_groups.copy()
    targets_comp.append('Water_Ions')
    index_ori = (gmxdir / "index.ndx")
    if index_ori.is_file():
        logger.debug(f'Index file {gmxdir}/index.ndx already exists. Checking for target groups {targets_comp}')
        with index_ori.open('r') as f:
            ndx_groups = []
            for line in f:
                line = line.strip()
                if line.startswith("[") and line.endswith("]"):
                    ndx_groups.append(line[1:-1].strip())
            logger.debug(f'Groups found in existing index file: {ndx_groups}')
        if set(targets_comp).issubset(ndx_groups):
            logger.warning(f'The index file {gmxdir}/index.ndx already exists and contains the target groups {targets_comp}. Skipping index creation step.')
            return True
        else:
            logger.debug(f'The index file {gmxdir}/index.ndx does not contain all target groups {targets_comp}. Proceeding with index creation step.')
            return False
        
def create_index_file(logger, gmxdir: Path, structurefile: Path, dna: bool, rna: bool, phos: bool, mtl: list[str], dna_nt: list[str], rna_nt: list[str], chains=None, outfile: str = "index.ndx"):
    """
    Create an index file with the target groups.
    
    Arguments:
        logger        : Logger object.
        gmxdir        : Path to the directory with the MD results (e.g. rep1/gromacs).
        structurefile : Path to PDB structure file.
        dna           : True if DNA is present in the system.
        rna           : True if RNA is present in the system.
        phos          : True if the system contains phosphorylated residues.
        mtl           : List of unique metals present in the system.
        dna_nt        : List of unique DNA residue names present in the system.
        rna_nt        : List of unique RNA residue names present in the system.
        chains        : List of chain IDs for which to create groups (e.g. ['A', 'B'])
        outfile       : Path to output .ndx file.
    """

    PHOSRES = ["PTR", "TPO", "SEP"]

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

    # Extract group names
    group_pattern = re.compile(r"^\s*\d+\s+([A-Za-z0-9_+\-]+)\s*:", re.MULTILINE)
    groups = group_pattern.findall(index_output)
    group_n = len(groups)
    logger.debug(f"Found {group_n} groups in index file {gmxdir}/temp.ndx: {groups}")

    group_protein = groups.index("Protein")
    group_water = groups.index("Water")
    group_na = groups.index("Na")
    group_cl = groups.index("Cl")
    group_phosphate = []
    if dna and "DNA" in groups:
        group_dna = groups.index("DNA")
    if rna and "RNA" in groups:
        group_rna = groups.index("RNA")
    if phos:
        group_phosphate = [str(groups.index(group)) for group in groups if group in PHOSRES]
        logger.debug(f"list of phosphate groups found: {group_phosphate}")

    # Remove temp index file
    os.remove(index_temp)

    # Run again to create index file with target groups
    groups_input = ""
    i = group_n
    # Define groups
    if len(mtl) > 0:
        groups_input += f"""
        {"r " + " | r ".join(mtl)}
        name {i} Metals
        """
        group_metals = i
        i += 1
        groups_input += f"""
        {group_protein} | {group_metals}
        name {i} Protein_Metals
        """
        i += 1

    # make_ndx makes a group for only the first of the DNA/RNA/Other residue types, and does
    # not know some terminal names (e.g. AMBER A5/U3), so missing groups are built from names.
    for kind, wanted, names in (("DNA", dna, dna_nt), ("RNA", rna, rna_nt)):
        if not wanted or kind in groups:
            continue
        if not names:
            logger.error(f"{kind} is expected in the system, but make_ndx made no {kind} group and "
                         f"no known {kind} residue name was found in {structurefile}.")
            sys.exit(1)
        groups_input += f"""
        {"r " + " | r ".join(names)}
        name {i} {kind}
        """
        logger.debug(f"{kind} group built from residue names {names} as number {i}")
        if kind == "DNA":
            group_dna = i
        else:
            group_rna = i
        i += 1

    if phos and len(group_phosphate) > 0:
            groups_input += f"""
            {group_protein} | {" | ".join(group_phosphate)}
            name {i} Protein_P
            """
            i += 1

    # Define conditional group combinations
    parts = [group_protein]
    labels = ["Protein"]
    if dna:
        parts.append(group_dna)
        labels.append("DNA")
    if rna:
        parts.append(group_rna)
        labels.append("RNA")
    if phos:
        parts.extend(group_phosphate)
        labels.append("P")
    if len(mtl) > 0:
        parts.append(group_metals)
        labels.append("Metals")

    if len(parts) > 1:
        groups_input += f"""
        {' | '.join(str(part) for part in parts)}
        name {i} {'_'.join(labels)}
        """
        i += 1

    # Define chain groups if chains are specified
    logger.debug(f"Chains specified for group creation: {chains}")
    if isinstance(chains, list):
        for c in chains:
            if isinstance(c, str):
                logger.debug(f"Creating group for chain {c} as number {i}")
                groups_input += f"""chain {c}
                name {i} Chain_{c}
                """
                i += 1

    # Water and ions (always present, + quit at the end)
    groups_input += f"""
        {group_water} | {group_na} | {group_cl}
        name {i} Water_Ions
        quit
        """
    i += 1

    logger.debug(f"Final group input: {groups_input}")

    cmd = [ "gmx", "make_ndx", 
        "-f", structurefile,
        "-o", f"{gmxdir}/temp_groups.ndx",
        "-quiet"
        ]
    process = subprocess.run(cmd, input=groups_input, text=True, capture_output=True, check=True)
    log_process_output(logger, process, "gmx make_ndx")
    logger.info("Index file created with new groups")
    # Move temp_groups.ndx to final index.ndx
    os.replace(f"{gmxdir}/temp_groups.ndx", f"{gmxdir}/{outfile}")
    if (gmxdir / f"{outfile}").is_file():
        logger.info(f"Index file with target groups saved as {gmxdir}/{outfile}")
    return 

def group_index(gmxdir: Path = None, logfile: Path = None, logname: str = 'log', structurefile: Path = None, dna: bool = False, rna: bool = False, phos: bool = False):
    """
    Get the group name of the main complex, insert it into .mdp parameter files 
    and call create_index_file to generate index with required groups.
    
    Arguments:
        gmxdir        : Path to the directory with the MD results (e.g. rep1/gromacs).
        logfile   : Path to the log file.
        logname   : Name of the logger.
        structurefile : Path to PDB structure file.
        dna           : True if DNA is present in the system.
        rna           : True if RNA is present in the system.
        phos          : True if the system contains phosphorylated residues.
    """    
    ##### INITIALIZE LOG
    logger = setup_logger(logname, logfile)

    # Get main group name and list of metals and nucleotides present in the structure file
    group_main, group_protein, mtl, dna_nt, rna_nt = get_name_main(logger, structurefile, dna, rna, phos)

    # In the MD step, substitute group name in nvt.mdp, npt.mdp, md.mdp
    if logname == 'md' and group_main != 'Protein':
        for mdp_file in ['nvt.mdp', 'npt.mdp', 'md.mdp']:
            mdp_path = gmxdir / mdp_file
            with mdp_path.open('r') as f:
                contents = f.read()
            contents = re.sub(r'\S*Protein\S*', group_main, contents)
            with mdp_path.open('w') as f:
                f.write(contents)

    # Check if an index file with the target groups already exists
    if targets_in_index(logger, gmxdir, [group_main]):
        return group_main
    else: # If not, create index file with target groups
        create_index_file(logger, gmxdir, structurefile, dna, rna, phos, mtl, dna_nt, rna_nt)
    
    return group_main


def group_index_interactions(gmxdir: str = None, logfile: str = None, logname: str = 'log', structurefile: str = None, indexfile: str = None,
                             group1: str = None, group2: str = None, dna: bool = False, rna: bool = False, phos: bool = False, chains: dict = None): 
    """
    Get the group name of the main complex, insert it into .mdp parameter files 
    and call create_index_file to generate index with required groups.
    
    Arguments:
        gmxdir        : Path to the directory with the MD results (e.g. rep1/gromacs).
        logfile   : Path to the log file.
        logname   : Name of the logger.
        structurefile : Path to PDB structure file.
        indexfile : Path to output .ndx file to generate, e.g. interaction_group1_group2.ndx
        group1    : First atom group for hydrogen bond analysis.
        group2    : Second atom group for hydrogen bond analysis.
        dna           : True if DNA is present in the system.
        rna           : True if RNA is present in the system.
        phos          : True if the system contains phosphorylated residues.
        chains          : Dictionary of original ID chains and corresponding updated IDs, and offsets.
    """    
    
    ##### INITIALIZE LOG
    logger = setup_logger(logname, logfile)

    # Check that group1 and group2 are provided
    if group1 is None or group2 is None:
        logger.error("Both group1 and group2 must be specified for interaction analysis.")
        sys.exit(1)

    # Get main group name and list of metals and nucleotides present in the structure file
    group_main, group_protein, mtl, dna_nt, rna_nt = get_name_main(logger, structurefile, dna, rna, phos)

    # Check if an index file with the target groups already exists
    ndx_names = []
    target_chains = []
    for i, group in enumerate([group1, group2]):
            if re.match(r"^[A-Za-z]$", group):
                # Map old chain IDs to new ones
                new_id = next((k for k, v in chains.items() if v['id_old'] == group), None)
                if new_id:
                    ndx_names.append(f"Chain_{new_id}")
                    target_chains.append(new_id)
            elif group == "DNA" and dna:
                ndx_names.append("DNA")
            elif group == "RNA" and rna:
                ndx_names.append("RNA")
            elif group == "Protein":
                ndx_names.append(group_protein)
    logger.debug(f"Target groups for interaction analysis: {ndx_names}")
    # Check that ndx_names contains two names
    if len(ndx_names) != 2:
        logger.error(f"Error in creating list of index group names. Expected 2 names, got {len(ndx_names)}.")
        sys.exit(1)

    
    # If index file with target groups exists, create copy with target group names in the title 
    if targets_in_index(logger, gmxdir, ndx_names):
        shutil.copy(gmxdir / "index.ndx", gmxdir / f"interaction_{group1}_{group2}.ndx")
        logger.debug(f"Copied index as {gmxdir / f'interaction_{group1}_{group2}.ndx'}")
    # If not, create index file with target groups
    else: 
        create_index_file(logger, gmxdir, structurefile, dna, rna, phos, mtl, dna_nt, rna_nt, target_chains, f"interaction_{group1}_{group2}.ndx")
    logger.debug(f"Returning {ndx_names}")
    return ndx_names
