#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import glob
import time
from scripts.utils.helpers import setup_logger
from pathlib import Path

def main(logfile: Path, gmxdir: Path, topol: Path = Path("./topol.top"), range_a: list = [0, -1], range_b: list = [0, -1],
         type: int = 1, type2: int = 2, low: float = 0.0, up1: float = 0.3, up2: float = 0.4, fac: float = 1.0):
    """
    Calculate distance restraints for DNA H-bonds and append them to a topology file.

    Arguments:
        logfile : Path to the log file.
        gmxdir  : Directory containing GROMACS files
        topol   : Topology file to append restraints to (default: ./topol.top)
        range_a : Range of residues for DNA A (0-indexed) (default: [0, -1], full length)
        range_b : Range of residues for DNA B (0-indexed) (default: [0, -1], full length)
        type    : Restraint type (default: 1)
        type2   : Secondary restraint type (default: 2)
        low     : Lower bound distance (Å) (default: 0.0)
        up1     : Upper bound distance 1 (Å) (default: 0.3)
        up2     : Upper bound distance 2 (Å) (default: 0.4)
        fac     : Force constant scaling factor (default: 1.0)
    """

    ##### Initialize logger
    logger = setup_logger("prep", logfile)
    toppar = gmxdir / "toppar"

    ##### CHECK INPUT #####
    if not gmxdir.is_dir():
        logger.error(f"Directory {str(gmxdir)} not found")
        sys.exit(1)
    if not toppar.is_dir():
        logger.error(f"Topology directory {str(toppar)} not found")
        sys.exit(1)
    if not topol.is_file():
        logger.error(f"Topology file {str(topol)} not found.")
        sys.exit(1)

    # Check if checkpoint file from previous run exists
    checkpoint_files = list(gmxdir.glob(f"check_05_restraindnahbonds_*.done"))
    if len(checkpoint_files) > 0:
        logger.warning("Checkpoint file found.")
        # Check if topol file contains distance_restraints section
        with open(topol, 'r') as f:
            if "[ distance_restraints ]" in f.read():
                logger.info(f"[ distance_restraints ] section found in {str(topol)}.")
                logger.warning("Skipping DNA H-bond restraints addition.")
                return
            else:
                logger.warning(f"Checkpoint file found but [ distance_restraints ] section not found in {str(topol)}. The previous run might have been incomplete. Proceeding with DNA H-bond restraints addition.")

    # Find DNA itp files
    dna_resnames = {"DC", "DA", "DG", "DT"}
    dna_files = []
    # Find DNA itp files
    for file in toppar.glob("*.itp"):
        has_dna = False
        with open(file) as f:
            for line in f:
                line = line.strip()
                # match comment lines that define residues
                # example line:
                # ; residue    2 DG rtp DG q -1.0
                if line.startswith(";") and "residue" in line:
                    cols = line.split()
                    if len(cols) >= 3:
                        resname = cols[3]
                        if resname in dna_resnames:
                            has_dna = True
                            dna_files.append(Path(file))
                            logger.debug(f"Found DNA residue {resname} in {file}.")
                            break
    if len(dna_files) != 2:
        logger.error(f"--dna flag provided but found {len(dna_files)} files containing DNA residues. Expected 2 files.")
        sys.exit(1)
    else:
        dna_a, dna_b = dna_files
        logger.info("The complex includes DNA. Adding DNA H-bond restraints.")
    
    if not dna_a.is_file():
        logger.error(f"DNA chain A topology file {str(dna_a)} not found.")
        sys.exit(1)
    if not dna_b.is_file():
        logger.error(f"DNA chain B topology file {str(dna_b)} not found.")
        sys.exit(1)

    # distance bounds order
    if low > up1 or up1 > up2:
        logger.error("Invalid distance bounds: ensure that low <= up1 <= up2.")
        sys.exit(1)

    # fac should be positive
    if fac <= 0:
        logger.error("Force constant scaling factor (fac) should be positive.")
        sys.exit(1)

    if range_a is None:
        range_a = [0, -1]
    if range_b is None:
        range_b = [0, -1]

    ##### CHECK DNA SEQUENCES #####
    # Extract lines between [atoms] and [bonds]
    seqs = {}
    dna_files = {'a': dna_a, 'b': dna_b}
    dna_dicts = {'a': [], 'b': []}
    for key, dna_file in dna_files.items():
        with open(dna_file) as f:
            dna_atoms_lines = []
            resnr = []
            inside = False
            for line in f:
                if "[ atoms ]" in line:
                    inside = True
                    continue
                if "[ bonds ]" in line:
                    inside = False
                if inside and not line.lstrip().startswith(";"):
                    dna_atoms_lines.append(line)
            seq = []
            for line in dna_atoms_lines:
                cols = line.split()
                if len(cols) < 5:
                    continue
                # Extract sequence
                if cols[2] not in resnr:
                    resnr.append(cols[2])
                    seq.append(cols[3][1]) # substr($4,2,1) → second character of parts[3]
                # Save columns resnr, residue, atom to dict
                select_cols = cols[2:6] 
                dna_dicts[key].append(select_cols)
        seqs[key] = seq

    for key, value in seqs.items():
        logger.debug(f"Sequence for DNA {key.upper()}: {''.join(seqs[key])}")

    # Apply residue range filters if provided
    ranges  = [range_a, range_b]
    for i, r in enumerate(ranges):
        if int(r[0]) != 0 or int(r[1]) != -1:
            key = 'a' if i == 0 else 'b'
            logger.debug(f"Applying residue range filter for DNA {key.upper()}: {r}")
            seqs[key] = seqs[key][int(r[0]) : int(r[1]) if int(r[1]) != -1 else None]
            logger.debug(f"Filtered sequence for DNA {key.upper()}: {''.join(seqs[key])}")

    # Check equal lengths (in the selected ranges)
    len_a = len(seqs['a'])
    len_b = len(seqs['b'])
    if len_a == len_b:
        logger.info(f"DNA sequences have the same length {len_a}.")
    else:
        logger.error(f"DNA sequences have different lengths: {len_a} and {len_b}.")
        sys.exit(1)

    # Check valid DNA sequences
    nts = {'A', 'T', 'C', 'G'}
    for key, seq in seqs.items():
        for nt in seq:
            if nt not in nts:
                logger.error(f"Invalid nucleotide '{nt}' found in DNA {key.upper()} sequence. Only A, C, G, T are allowed.")
                sys.exit(1)
    logger.info("Both DNA sequences are valid.")

    # Check complementary seqs
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for i in range(len_a):
        if seqs['a'][i] != complement[seqs['b'][len_b - i - 1]]:
            logger.error(f"Nucleotides at position {i+1} are not complementary: {seqs['a'][i]}-{seqs['b'][len_b - i - 1]}.")
            sys.exit(1)
    logger.info("DNA sequences are complementary.")

    ##### CALCULATE DISTANCE RESTRAINTS #####
    logger.debug(f"""Calculating distance restraints for h-bonds in {dna_a} and {dna_b}.\n
        Restrains will be writen to {topol}.
        Parameters: 
        type    index    type'    low    up1    up2    fac
        {type}	i	 {type2}	  {low}	 {up1}	{up2}   {fac}""")
    
    with open(topol, 'a') as f:
        text = """
[ intermolecular_interactions ]
[ distance_restraints ]
; ai    aj    type    index    type'    low    up1    up2    fac\n"""
        f.write(text)

    # Mapping local ids to global ids
    with open(topol) as f:
        itps = []
        molmap = {}
        in_molecules = False
        for line in f:
            # Get list of itps included in topol (except forcefield)
            if line.strip().startswith("#include") and "forcefield" not in line:
                filename = line.split("/")[-1].strip().strip('"') 
                itps.append(filename)
            # Get number of copies per moleculetype
            if line.strip().startswith('[ molecules ]'):
                in_molecules = True
                continue
            # End of [ molecules ] section
            if line.strip().startswith("[ intermolecular_interactions ]"):
                in_molecules = False
                continue
            if in_molecules:
                if not line or len(line) == 0 or line.startswith(';'):
                    continue
                # First data line, i.e. "system1          2"
                cols = line.split()
                if len(cols) == 2:
                    molmap[cols[0]] = int(cols[1])
                elif len(cols) > 0:
                    logger.warning(f"[ molecules ] line {line} in {topol} could not be split in columns. Skipping count of this compound.")


    global_id = 0
    global_a = 0
    global_b = 0
    logger.debug(f"Found .itp files included in {topol} (excluding forcefield.itp):")


    for itp in itps:
        logger.debug(itp)
        if dna_a.name == itp:
            logger.debug(f"Match: {itp}")
            global_a = global_id
        elif dna_b.name == itp:
            logger.debug(f"Match: {itp}")
            global_b = global_id
        in_atoms = False
        in_moleculetype = False
        itp_name = None
        with open(f"{gmxdir}/toppar/{itp}") as f:
            last_value = None
            for line_raw in f:
                line = line_raw.strip()
                # ---- [ moleculetype ] handling ----
                if line.startswith('[') and 'moleculetype' in line:
                    in_moleculetype = True
                    continue
                if in_moleculetype and itp_name is None:
                    if not line or line.startswith(';'):
                        continue
                    # First data line, i.e. "system1          3"
                    cols = line.split()
                    if cols:
                        itp_name = cols[0]
                    in_moleculetype = False
                # Find [ atoms ] section
                if "[ atoms ]" in line:
                    in_atoms = True
                    continue
                if "[ bonds ]" in line:
                    in_atoms = False
                if not in_atoms:
                    continue
                cols = line.split()
                if cols and cols[0].isnumeric():
                        last_value = cols[0]
            logger.debug(f"File {itp} corresponds to compound {itp_name}, containing {last_value} atoms.")
            if itp_name and itp_name in molmap.keys():
                natoms = int(last_value) * molmap[itp_name]
                logger.debug(f"Found {molmap[itp_name]} ocurrences of compound {itp_name} in the system.")
            else:
                natoms = int(last_value)
                logger.debug(f"Number of ocurrences of compound {itp_name} in the system could not be determined. Using default 1.")
        # after reading whole itp file:
        if natoms is not None:
            global_id += natoms
            

    logger.debug(f"Number of atoms in the complex before DNA A: {global_a}")
    logger.debug(f"Number of atoms in the complex before DNA B: {global_b}")

    # Map atom pairs by name
    pairs = {'A': {'atoms_a':['N6', 'N1'], 'atoms_b': ['O4', 'N3'], 'comp': 'DT'},
             'T': {'atoms_a':["O4", "N3"], 'atoms_b':["N6", "N1"], 'comp': 'DA'},
             'C': {'atoms_a':["O2", "N3", "N4"], 'atoms_b':["N2", "N1", "O6"], 'comp': 'DG'},
             'G': {'atoms_a':["N2", "N1", "O6"], 'atoms_b':["O2", "N3", "N4"], 'comp': 'DC'}}
    for i, nt in enumerate(seqs['a']):
        nt_pairs = pairs[nt]
        atoms_a = nt_pairs['atoms_a']
        atoms_b = nt_pairs['atoms_b']
        comp = nt_pairs['comp']
        with open(topol, "a") as out:
            for idx, atom in enumerate(atoms_a):
                i_val = i + 1 + int(range_a[0])
                comp_i = len_a - i + int(range_b[0])

                for row in dna_dicts['a']:
                    if row[0] == str(i_val) and f"D{nt}" in row[1] and row[2] == atoms_a[idx]:
                        atom_a = int(row[3])
                        break
                for row in dna_dicts['b']:
                    if row[0] == str(comp_i) and comp in row[1] and row[2] == atoms_b[idx]:
                        atom_b = int(row[3])
                        break
                # Add atom pair to distance restrains table
                out.write(
                    f"{atom_a + global_a}    {atom_b + global_b}    {type}    {i + 1}    {type2}    {low}    {up1}    {up2}    {fac}\n")

    logger.info(f"Distance restraints for DNA H-bonds added to {topol}.")
    # Create checkpoint file
    checkpoint_filename = gmxdir / f"check_05_restraindnahbonds_{int(time.time())}.done"
    with open(checkpoint_filename, 'w') as f:
        f.write("Distance restraints for DNA H-bonds added successfully.\n")


if __name__ == "__main__":
    main()