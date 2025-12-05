#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import time
import sys
import shutil
import glob
import argparse
from scripts.utils.helpers import setup_logger


def remove_hydrogens(dir, filename, logger_name, logfile):
	logger = setup_logger(logger_name, logfile)
	cmd = f"""
	pdb4amber -i {dir}/{filename}.pdb -o {dir}/{filename}_noH.pdb --nohyd --no-conect
	pdb4amber -i {dir}/{filename}_noH.pdb -o {dir}/{filename}_reduced.pdb --reduce --no-conect
    """
	process = subprocess.run(cmd, shell=True, text=True, capture_output=True)
	if process.returncode != 0:
		logger.error(process.stderr)
		return None  
	else:
		logger.info(f"Hydrogen processing complete.\n {process.stdout}")
	return f"{dir}/{filename}_reduced.pdb"


def add_hydrogens(file, logger_name, logfile, ffs=['water.opc', 'protein.ff19SB']):
	"""Add hydrogen atoms to a PDB file using tleap."""
	logger = setup_logger(logger_name, logfile)
	protonated_file = file.replace('_reduced.pdb', '_protonated.pdb')
	# Build tleap input dynamically
	ffs_import = "\n".join(f"source leaprc.{ff}" for ff in ffs)

	tleap_input = f"""
		{ffs_import}
		complex = loadPdb {file}
		savePdb complex {protonated_file}
		quit
	"""
	logger.debug(f"Running tleap with input:\n{tleap_input}")  # Debugging: Show tleap input
	process = subprocess.run(['tleap', '-f', '-'], input=tleap_input, text=True, capture_output=True)
	
	if process.returncode != 0:
		logger.error(f"tleap error:\n {process.stderr}")
	else:
		logger.debug(f"tleap output:\n {process.stdout}")

	return protonated_file


def gen_chains(complex_file, logger_name, logfile, ffs=['water.opc', 'protein.ff19SB'], *chains):
	"""Generate separate PDB files.
		*chains takes a varible number of arguments for the chains to be separated, for example:
		mono_a='A', mono_b='B', dna_a='C', dna_b='D'
	"""
	logger = setup_logger(logger_name, logfile)
	complex_name = os.path.basename(complex_file)

	files = []

	for chain in chains:
		chain_file = complex_file.replace(complex_name, f'{chain}_noH.pdb')

		# Run pdb_selchain.py command for chain
		cmd_chain = ['pdb_selchain', f'-{chain}', complex_file]
		process_chain = subprocess.run(cmd_chain, capture_output=True, text=True)

		# Check for errors
		if process_chain.returncode != 0:
			logger.error(f"Error occurred while generating receptor file:\n {process_chain.stderr}")
			return None
		else:
			logger.info(f"Chain {chain} extraction complete.")

		# Write receptor data to file
		with open(chain_file, 'w') as f:
			f.write(process_chain.stdout)

		chain_file = add_hydrogens(chain_file, logger_name, ffs)
		files.append(chain_file)

	return tuple(files)


def solvate_gen_parm(complex_file, ffs=['water.opc', 'protein.ff19SB'], *chains):
	"""Protonate and solvate the complex and chains using LEaP."""
	complex_name = os.path.splitext(complex_file)[0]
	complex_solvated_file = complex_name + '_solvated.pdb'

	filebases = []
	for chain in chains:
		filebases.append(chain[:-4])


	# Build tleap input dynamically
	ffs_import = "\n".join(f"source leaprc.{ff}" for ff in ffs)

	cmd = f'''
	{ffs_import}
	loadamberparams frcmod.ionslm_1264_opc

	complex = loadPdb {complex_file}
	'''

	# Dynamically add loadPdb commands
	for chain in filebases:
		cmd += f'{chain} = loadPdb {chain}.pdb\n'

	cmd += '\nset default PBRadii mbondi2\n\n'

	# Save Amber prmtop/inpcrd for complex
	cmd += f'saveAmberParm complex {complex_name}.prmtop {complex_name}.inpcrd\n'

	# Dynamically save Amber prmtop/inpcrd for all mono_ and dna_ files
	for chain in filebases:
		cmd += f'saveAmberParm {chain} {chain}.prmtop {chain}.inpcrd\n'

	# Solvate and add ions
	cmd += f'''
	solvateOct complex OPCBOX 10.0
	addIons2 complex Na+ 1  Cl- 1
	saveAmberParm complex {complex_name}_solvated.prmtop {complex_name}_solvated.inpcrd
	savePdb complex {complex_solvated_file}

	quit
	'''
	# Encode the command string into bytes
	cmd_bytes = cmd.encode('utf-8')

	# Run the command using subprocess
	subprocess.run(['tleap', '-f', '-'], input=cmd_bytes)




########## MAIN ##########

def main():
	start_time = time.time()

	parser = argparse.ArgumentParser(description="Prepare system using AmberTools. Add solvent and ions, generate prmtop and inpcrd files based on force fields.")
	parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
	parser.add_argument('--path-pdb', dest='path_pdb', help='Input pdb file path', required=True)
	parser.add_argument('--dna', dest='with_dna', action='store_true', help='Use if the system contains DNA. Adds DNA.bsc1 force field.')
	parser.add_argument('--phos', dest='with_phosphate', action='store_true', help='Use if the system contains phosphorylated residues. Adds phosaa19SB force field.')
	parser.add_argument('--ffs', dest='ffs', help='Comma separated list of force fields to use from the list available in LEaP. Avoid using clashing force fields. Default: protein.ff19SB,water.opc')

	# print help if no arguments are provided
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(0)

	args = parser.parse_args()

	amberdir, pdb = os.path.split(args.path_pdb)
	filename = os.path.splitext(pdb)[0]

	ff_default = ['water.opc']

	##### Initialize logger
	logger_name = 'prep'
	logger = setup_logger(logger_name, args.logfile)

	##### Check for checkpoint files from previous runs #####
	expected_files = [
		f"{amberdir}/{filename}_protonated_solvated.inpcrd",
		f"{amberdir}/{filename}_protonated_solvated.prmtop",
	]

	checkpoint_files = glob.glob(f"{amberdir}/check_01_amber_*.done")
	if checkpoint_files:
		logger.warning(f"Checkpoint file(s) found: {', '.join(checkpoint_files)}")
		# Check if expected output files exist
		missing = [p for p in expected_files if not os.path.isfile(p)]
		if not missing:
			logger.warning("The expected output files are already present. Skipping AmberTools preparation.")
			sys.exit(0)


	# Check if input file is a PDB file
	if not pdb.lower().endswith('.pdb'):
		logger.error("Input file must be a PDB file with .pdb extension.")
		sys.exit(1)

	if args.with_dna:
		logger.debug("Flag --dna used. This is a DNA-protein complex.")
		ff_default.append('DNA.bsc1')

	if args.with_phosphate:
		logger.debug("Flag --phos used. This system contains phosphorylated residues.")
		ff_default.append('phosaa19SB')
	else:
		ff_default.append('protein.ff19SB')

	if args.ffs:
		logger.debug(f"Custom force fields provided: {args.ffs}")
		ffs = args.ffs.split(",")
	else:
		logger.debug(f"No custom force fields provided. Using default: {ff_default}")
		ffs = ff_default

	# Prepare the complex (make sure the hydrogens naming matches amber)
	noH_complex = remove_hydrogens(amberdir, filename, logger_name, args.logfile)
	protonated_complex = add_hydrogens(noH_complex, logger_name, args.logfile, ffs)

	# Generate the receptor & ligand (also match the hydrogens to amber specifics)
	if args.with_dna:
		mono_a_file, mono_b_file, dna_a_file, dna_b_file = gen_chains(noH_complex, logger_name, args.logfile, ffs,  'A', 'B', 'C', 'D')
		solvate_gen_parm(protonated_complex, ffs, mono_a_file, mono_b_file, dna_a_file, dna_b_file)
	else:
		mono_a_file, mono_b_file = gen_chains(noH_complex, logger_name, args.logfile, ffs, 'A', 'B')	
		solvate_gen_parm(protonated_complex, ffs, mono_a_file, mono_b_file)

	solvated_complex_prmtop = protonated_complex.replace('.pdb', '_solvated.prmtop')
	solvated_complex_pdb = protonated_complex.replace('.pdb', '_solvated.pdb')

	end_time = time.time()
	elapsed_time = end_time - start_time
	logger.debug(f"Total execution time: {elapsed_time:.2f} seconds")

	# Move logs to amberdir
	if os.path.isfile("leap.log"):
		shutil.move("leap.log", f"{amberdir}/leap.log")
	if os.path.isfile("reduce_info.log"):
		shutil.move("reduce_info.log", f"{amberdir}/reduce_info.log")

	# Check output files
	logger.debug("Checking generated files.")

	missing = [p for p in expected_files if not os.path.isfile(p)]
	for file in expected_files:
		if file not in missing:
			logger.info(f"generated {file}")
		else:
			logger.error(f"{file} not generated")

	if not missing:
		logger.info("AmberTools preparation completed successfully.")
		checkpoint_filename = f"{amberdir}/check_01_amber_{int(time.time())}.done"
		with open(checkpoint_filename, 'w') as f:
			f.write("AmberTools preparation completed successfully.\n")
		sys.exit(0)
	else:
		logger.error("AmberTools preparation failed. There are missing files.")
		sys.exit(1)


if __name__ == "__main__":
	main()


