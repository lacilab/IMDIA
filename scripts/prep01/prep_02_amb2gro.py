#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import glob
import time
import parmed as pmd
import argparse
from scripts.utils.colors import ok, warn, error

def main():
    parser = argparse.ArgumentParser(description="Convert Amber files to GROMACS format.")
    parser.add_argument("-p", "--prmtop", dest="prmtop", type=str, required=True,
                        help="Prmtop file")
    parser.add_argument("-c", "--inpcrd", dest="inpcrd", type=str, required=True,
                        help="Inpcrd file")
    parser.add_argument("-t", "--top", dest="top", type=str, required=True,
                        help="GROMACS top file")
    parser.add_argument("-g", "--gro", dest="gro", type=str, required=True,
                        help="GROMACS gro file")

    # print help if no arguments are provided
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    amberdir = os.path.dirname(os.path.abspath(args.top))

    ##### Check for checkpoint files from previous runs #####
    expected_files = [
        f"{args.top}",
        f"{args.gro}"
    ]

    checkpoint_files = glob.glob(f"{amberdir}/check_02_amb2gro_*.done")
    if checkpoint_files:
        warn(f"Checkpoint file(s) found: {', '.join(checkpoint_files)}")
        # Check if expected output files exist
        missing = [p for p in expected_files if not os.path.isfile(p)]
        if not missing:
            warn("The expected output files are already present. Skipping Amber to GROMACS conversion.")
            sys.exit(0)

    # If checkpoint files are not found, proceed with the conversion
    amber = pmd.load_file(args.prmtop,args.inpcrd)

    # Cmap correction, including PTR
    CX_dir = { 'ALA' : 'XC0', 'ARG' : 'XC1', 'ASH' : 'XC2', 'ASN' : 'XC3', 'ASP' : 'XC4', 'CYM' : 'XC5', 'CYS' : 'XC6', 
            'CYX' : 'XC5', 'GLH' : 'XC7', 'GLN' : 'XC8', 'GLU' : 'XC9', 'GLY' : 'XC10', 'HID' : 'XC5', 'HIE' : 'XC5',
            'HIP' : 'XC5', 'HYP' : 'XC5', 'ILE' : 'XC11', 'LEU' : 'XC5', 'LYN' : 'XC12', 'LYS' : 'XC12', 'MET' : 'XC5',
            'PHE' : 'XC5', 'PRO' : 'XC13', 'SER' : 'XC14', 'THR' : 'XC15', 'TRP' : 'XC5', 'TYR' : 'XC5', 'VAL' : 'XC11',
                'PTR' : 'XC5' } 

    known_types = {}
    known_types_list = []
    for cmap in amber.cmaps:  # loop through all cmap terms
        atom_type = cmap.atom3.atom_type # get atom type
        new_typename = CX_dir[cmap.atom3.residue.name] # get new CA type name (e.g. XC0)
        if new_typename in known_types_list: # if new type name in known types, directly assign atom to new type
            cmap.atom3.atom_type = known_types[cmap.type]
            cmap.atom3.type = known_types[cmap.type].name
            continue # skip

        # if new type name is not in known types, create new type
        new_type = pmd.AtomType(new_typename, atom_type.number, atom_type.mass, atom_type.atomic_number, new_typename, atom_type.charge)
        new_type.epsilon = atom_type.epsilon # copy over epsilon and rmin
        new_type.rmin = atom_type.rmin
        known_types[cmap.type] = new_type
        known_types_list.append(new_typename)

        # assgin atom to new type
        cmap.atom3.atom_type = known_types[cmap.type]
        cmap.atom3.type = known_types[cmap.type].name


    # output 
    amber.save(args.top, overwrite=True, format='gromacs')
    amber.save(args.gro, overwrite=True, format='gro')

    print("\nChecking generated files.")

    missing = [p for p in expected_files if not os.path.isfile(p)]
    for file in expected_files:
        if file not in missing:
            ok(f"generated {file}")
        else:
            error(f"{file} not generated")

    if not missing:
        ok("Amber to GROMACS conversion completed successfully.")
        checkpoint_filename = f"{amberdir}/check_02_amb2gro_{int(time.time())}.done"
        with open(checkpoint_filename, 'w') as f:
            f.write("Amber to GROMACS conversion completed successfully.\n")
        sys.exit(0)
    else:
        error("Amber to GROMACS conversion failed. There are missing files.")
        sys.exit(1)

if __name__ == "__main__":
    main()