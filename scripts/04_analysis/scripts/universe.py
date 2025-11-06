#!/usr/bin/env python3

import sys
import MDAnalysis as mda
import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')
import os
import pickle
import argparse

# Arguments
parser = argparse.ArgumentParser(description="Create MDAnalysis Universe from PDB and XTC files.")
parser.add_argument("-f", "--pdb", help="Path to the PDB file")
parser.add_argument("-x", "--xtc", help="Path to the XTC file")
parser.add_argument("-d", "--outdir", help="Output directory")
parser.add_argument("-u", "--universe", help="Path to save the Universe object")

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(0)

args = parser.parse_args()
pdbfile = args.pdb
xtcfile = args.xtc
outdir = args.outdir
universe = args.universe

os.makedirs(outdir, exist_ok=True)

u = mda.Universe(pdbfile, xtcfile)
with open(universe, 'wb') as f:
    pickle.dump(u, f)

print(u)
print('Number of frames in the trajectory: ', len(u.trajectory))
print('Number of residues: ', len(u.residues))
print(u.residues)
print('Number of atoms: ', len(u.atoms))
print(u.atoms)
print('Segments: ', u.segments)

# Select atoms by chain
dimer = u.select_atoms("segid A B")  # Assuming dimer spans chains A and B
monoA = u.select_atoms("segid A")    # Monomer A is chain A
monoB = u.select_atoms("segid B")    # Monomer B is chain B
dna = u.select_atoms("segid C D")    # Assuming dsDNA spans chains C and D
dnaC = u.select_atoms("segid C")     # DNA chain A is chain C
dnaD = u.select_atoms("segid D")     # DNA chain B is chain D

# Store selections in a list for easy iteration
selections = [
    ("Dimer", dimer, "A, B"),
    ("Monomer A", monoA, "A"),
    ("Monomer B", monoB, "B"),
    ("dsDNA", dna, "C, D"),
    ("DNA A", dnaC, "C"),
    ("DNA B", dnaD, "D")
]
# Print atom group name, chain, and atom index range
for name, group, chain in selections:
    atom_indices = group.atoms.indices
    index_range = f"{min(atom_indices)}:{max(atom_indices)}"
    print(f"{name}: Chain(s) {chain}, Atom Index Range {index_range}, Residues {len(group.residues)}")