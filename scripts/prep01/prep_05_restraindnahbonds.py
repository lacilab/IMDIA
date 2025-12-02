#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import glob
import time
from scripts.utils.colors import ok, warn, error

def main():
    parser = argparse.ArgumentParser(description="Calculate distance restraints for DNA H-bonds and append them to a topology file.")
    parser.add_argument("--gmxdir", '-g', dest='gmxdir', type=str,
                        help='Directory containing GROMACS files')
    parser.add_argument("--topol", "-t", dest="topol", type=str, default="./topol.top",
                        help="Topology file to append restraints to")
    parser.add_argument("--dna_a", dest="dna_a", type=str, default="./toppar/molecule_3.itp",
                        help="DNA chain A topology file (.itp)")
    parser.add_argument("--dna_b", dest="dna_b", type=str, default="./toppar/molecule_4.itp",
                        help="DNA chain B topology file (.itp)")
    parser.add_argument("--range_a", dest="range_a", type=str, default="[0,-1]",
                        help="Range of residues to consider for DNA A. 1-indexed (e.g., [1,10]). Default: full length.")
    parser.add_argument("--range_b", dest="range_b", type=str, default="[0,-1]",
                        help="Range of residues to consider for DNA B. 1-indexed (e.g., [1,10]). Default: full length.")
    parser.add_argument("--type", dest="type", type=int, default=1,
                        help="Restraint type. Default: 1")
    parser.add_argument("--type2", dest="type2", type=int, default=2,
                        help="Secondary restraint type. Default: 2")
    parser.add_argument("--low", dest="low", type=float, default=0.0,
                        help="Lower bound distance (Å). Default: 0.0")
    parser.add_argument("--up1", dest="up1", type=float, default=0.3,
                        help="Upper bound distance 1 (Å). Default: 0.3")
    parser.add_argument("--up2", dest="up2", type=float, default=0.4,
                        help="Upper bound distance 2 (Å). Default: 0.4")
    parser.add_argument("--fac", dest="fac", type=float, default=1.0,
                        help="Force constant scaling factor. Default: 1.0")
    
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    ##### CHECK INPUT #####
    if not os.path.isdir(args.gmxdir):
        error(f"Directory {args.gmxdir} not found")
        sys.exit(1)
    if not os.path.isfile(args.topol):
        error(f"Topology file {args.topol} not found.")
        sys.exit(1)

    # Check if checkpoint file from previous run exists
    if glob.glob(f"{args.gmxdir}/check_05_restraindnahbonds_*.done"):
        warn("Checkpoint file found.")
        # Check if topol file contains distance_restraints section
        with open(args.topol, 'r') as f:
            if "[ distance_restraints ]" in f.read():
                ok(f"[ distance_restraints ] section found in {args.topol}.")
                warn("Skipping DNA H-bond restraints addition.")
                sys.exit(0)
            else:
                warn(f"Checkpoint file found but [ distance_restraints ] section not found in {args.topol}. The previous run might have been incomplete. Proceeding with DNA H-bond restraints addition.")

    if not os.path.isfile(args.dna_a):
        error(f"DNA chain A topology file {args.dna_a} not found.")
        sys.exit(1)
    if not os.path.isfile(args.dna_b):
        error(f"DNA chain B topology file {args.dna_b} not found.")
        sys.exit(1)

    # distance bounds order
    if args.low > args.up1 or args.up1 > args.up2:
        error("Invalid distance bounds: ensure that low <= up1 <= up2.")
        sys.exit(1)

    # fac should be positive
    if args.fac <= 0:
        error("Force constant scaling factor (fac) should be positive.")
        sys.exit(1)

    ##### CHECK DNA SEQUENCES #####
    # Extract lines between [atoms] and [bonds]
    seqs = {}
    dna_files = {'a': args.dna_a, 'b': args.dna_b}
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

    # Apply residue range filters if provided
    ranges  = {'a': args.range_a.split(','), 'b': args.range_b.split(',')}
    for key, value in ranges.items():
        if value[0] != 0 or value[1] != -1:
            print(f"Applying residue range filter for DNA {key.upper()}: {value}")
            seqs[key] = seqs[key][int(value[0]):int(value[1])]
            print(f"Filtered sequence for DNA {key.upper()}: {''.join(seqs[key])}")

    # Check equal lengths (in the selected ranges)
    len_a = len(seqs['a'])
    len_b = len(seqs['b'])
    if len_a == len_b:
        ok(f"DNA sequences have the same length {len_a}.")
    else:
        error(f"DNA sequences have different lengths: {len_a} and {len_b}.")
        sys.exit(1)

    # Check valid DNA sequences
    nts = {'A', 'T', 'C', 'G'}
    for key, seq in seqs.items():
        for nt in seq:
            if nt not in nts:
                error(f"Invalid nucleotide '{nt}' found in DNA {key.upper()} sequence. Only A, C, G, T are allowed.")
                sys.exit(1)
    ok("Both DNA sequences are valid.")

    # Check complementary seqs
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for i in range(len_a):
        if seqs['a'][i] != complement[seqs['b'][len_b - i - 1]]:
            error(f"Nucleotides at position {i+1} are not complementary: {seqs['a'][i]}-{seqs['b'][len_b - i - 1]}.")
            sys.exit(1)
    ok("DNA sequences are complementary.")

    ##### CALCULATE DISTANCE RESTRAINTS #####
    print(f"""Calculating distance restraints for h-bonds in {args.dna_a} and {args.dna_b}.\n
        Restrains will be writen to {args.topol}.
        Parameters: 
        type    index    type'    low    up1    up2    fac
        {args.type}	i	 {args.type2}	  {args.low}	 {args.up1}	{args.up2}   {args.fac}""")
    
    with open(args.topol, 'a') as f:
        text = """
[ intermolecular_interactions ]
[ distance_restraints ]
; ai    aj    type    index    type'    low    up1    up2    fac\n"""
        f.write(text)

    # Mapping local ids to global ids
    with open(args.topol) as f:
        itps = []
        for line in f:
            if line.startswith("#include") and "forcefield" not in line:
                filename = line.split("/")[-1].strip().strip('"') 
                itps.append(filename)
    
    global_id = 0
    global_a = 0
    global_b = 0
    print(f"Found .itp files included in {args.topol} (excluding forcefield.itp):")
    for itp in itps:
        print(itp)
        if os.path.basename(args.dna_a) in itp:
            print(f"Match: {args.dna_a}, {itp}")
            global_a = global_id
        elif os.path.basename(args.dna_b) in itp:
            print(f"Match: {args.dna_b}, {itp}")
            global_b = global_id
        n = None 
        inside = False
        with open(f"{args.gmxdir}/toppar/{itp}") as f:
            last_value = None
            for line in f:
                # Find [ atoms ] section
                if "[ atoms ]" in line:
                    inside = True
                    continue
                if "[ bonds ]" in line:
                    inside = False
                if not inside:
                    continue
                # Consider only lines with digits
                if not any(ch.isdigit() for ch in line):
                    continue
                # Split in columns
                parts = line.split()
                if not parts:
                    continue
                # take first column (atom id)
                last_value = parts[0]
        # after reading whole file:
        if last_value is not None:
            global_id += int(last_value)

    print(f"Number of atoms in the complex before DNA A: {global_a}")
    print(f"Number of atoms in the complex before DNA B: {global_b}")

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
        with open(args.topol, "a") as out:
            for idx, atom in enumerate(atoms_a):
                i_val = i + 1 + int(args.range_a.split(',')[0]) # 1-indexed
                comp_i = len_a - i

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
                    f"{atom_a + global_a}    {atom_b + global_b}    {args.type}    {i + 1}    {args.type2}    {args.low}    {args.up1}    {args.up2}    {args.fac}\n")

    ok(f"Distance restraints for DNA H-bonds added to {args.topol}.")
    # Create checkpoint file
    checkpoint_filename = f"{args.gmxdir}/check_05_restraindnahbonds_{int(time.time())}.done"
    with open(checkpoint_filename, 'w') as f:
        f.write("Distance restraints for DNA H-bonds added successfully.\n")


if __name__ == "__main__":
    main()