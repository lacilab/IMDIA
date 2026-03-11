#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import time
import sys
from pathlib import Path
from venv import logger
import MDAnalysis as mda
from scripts.utils.helpers import setup_logger


def remove_hydrogens(filename: str, logger_name: str, logfile: Path):
	logger = setup_logger(logger_name, logfile)
	logger.debug(f"Removing hydrogens from {filename} using pdb4amber.")
	outfile = filename.replace('.pdb', '_noH.pdb')

	cmd = f"""
	pdb4amber -i {filename} -o {outfile} --nohyd --no-conect
    """
	
	process = subprocess.run(cmd, shell=True, text=True, capture_output=True)
	if process.returncode != 0:
		logger.error(process.stderr)
		return None  
	else:
		logger.info(f"Hydrogen processing complete.\n {process.stdout}")
	return outfile


def add_hydrogens(file: str, logger_name: str, logfile: Path, ffs: list[str] = ['water.opc', 'protein.ff19SB']):
	"""Add hydrogen atoms to a PDB file using tleap."""
	logger = setup_logger(logger_name, logfile)
	logger.debug(f"Adding hydrogens to {file} using tleap with force fields: {ffs}")
	outfile = file.replace('_noH.pdb', '_protonated.pdb')
	# Build tleap input dynamically
	ffs_import = "\n".join(f"source leaprc.{ff}" for ff in ffs)

	tleap_input = f"""
		{ffs_import}
		complex = loadPdb {file}
		savePdb complex {outfile}
		quit
	"""
	logger.debug(f"Running tleap with input:\n{tleap_input}") 
	process = subprocess.run(['tleap', '-f', '-'], input=tleap_input, text=True, capture_output=True)
	
	if process.returncode != 0:
		logger.error(f"tleap error:\n {process.stderr}")
	else:
		logger.debug(f"tleap output:\n {process.stdout}")

	return outfile


def gen_chains(complex_file: Path, logger_name: str, logfile: Path, chains: list[str], ffs: list[str] = ['water.opc', 'protein.ff19SB']):
	"""
		Generate separate PDB files.
	"""
	logger = setup_logger(logger_name, logfile)
	logger.debug(f"Generating separate PDB files for chains: {chains} from {complex_file}")
	complex_name = Path(complex_file).stem

	files = []

	for chain in chains:
		# chain_file = complex_file.replace(complex_name, f'{chain}_noH.pdb')

		chain_file = complex_name + f'_{chain}_noH.pdb'

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

		chain_file = add_hydrogens(chain_file, logger_name, logfile, ffs)
		files.append(chain_file)

	return list(files)


def solvate_gen_parm(complex_file: Path, logger_name: str, logfile: Path, chains: list[str], ffs: list[str] = ['water.opc', 'protein.ff19SB']):
	"""Protonate and solvate the complex and chains using LEaP."""
	logger = setup_logger(logger_name, logfile)
	logger.debug(f"Solvating and generating Amber output for {complex_file} and chains: {chains}")
	complex_name = Path(complex_file).stem
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
	logger.debug(f"Running tleap with input:\n{cmd}")

	result = subprocess.run(['tleap', '-f', '-'], capture_output=True, text=True, check=True, input=cmd)
	lines = result.stdout.splitlines()
	filtered = []
	i = 0
	while i < len(lines):
		start = i
		while i < len(lines) and lines[i] == lines[start]:
			i += 1
		count = i - start
		if count > 10:
			filtered.append(f"SUPPRESSED {count} consecutive '{lines[start]}' lines.")
		else:
			filtered.extend(lines[start:i])
	logger.debug('\n'.join(filtered))
	




########## MAIN ##########

def main(logfile: Path = None, pdbpath: Path = None,
	dna : bool = False, phos : bool = False, ffs: list = None):
	start_time = time.time()

	# parser = argparse.ArgumentParser(description="Prepare system using AmberTools. Add solvent and ions, generate prmtop and inpcrd files based on force fields.")
	# parser.add_argument('--logfile', '-l', dest='logfile', default=None, help='Path to the log file.')
	# parser.add_argument('--path-pdb', dest='path_pdb', help='Input pdb file path', required=True)
	# parser.add_argument('--dna', dest='with_dna', action='store_true', help='Use if the system contains DNA. Adds DNA.bsc1 force field.')
	# parser.add_argument('--phos', dest='with_phosphate', action='store_true', help='Use if the system contains phosphorylated residues. Adds phosaa19SB force field.')
	# parser.add_argument('--ffs', dest='ffs', help='Comma separated list of force fields to use from the list available in LEaP. Avoid using clashing force fields. Default: protein.ff19SB,water.opc')
	
	##### Initialize logger
	logger_name = 'prep'
	logger = setup_logger(logger_name, logfile)

	amberdir = Path(pdbpath).parent
	pdbfile = Path(pdbpath).name
	structure = Path(pdbpath).stem

	logger.debug(f"Extracted amberdir: {amberdir}, pdbfile: {pdbfile}, structure: {structure}")


	ff_default = ['water.opc', 'protein.ff19SB']

	# Change to amberdir so fixed amber log files are stored here
	cwd = Path.cwd() 
	os.chdir(amberdir)


	##### Check for checkpoint files from previous runs #####
	expected_files = [
		Path(f"{structure}_protonated_solvated.inpcrd"),
		Path(f"{structure}_protonated_solvated.prmtop"),
	]

	checkpoint_files = list(Path.cwd().glob("check_01_amber_*.done"))
	if len(checkpoint_files) > 0:
		logger.warning(f"Checkpoint file(s) found: {', '.join(f.name for f in checkpoint_files)}")
		# Check if expected output files exist
		missing = [p for p in expected_files if not p.exists()]
		if not missing:
			logger.warning("The expected output files are already present. Skipping AmberTools preparation.")
			os.chdir(cwd)
			return

	# Check if input file is a PDB file
	if not pdbpath.suffix.lower() == '.pdb':
		logger.error("Input file must be a PDB file with .pdb extension.")
		sys.exit(1)

	if dna:
		logger.debug("Flag --dna used. This is a DNA-protein complex.")
		ff_default.append('DNA.bsc1')

	if phos:
		logger.debug("Flag --phos used. This system contains phosphorylated residues.")
		ff_default.append('phosaa19SB')
	#else:
	#	ff_default.append('protein.ff19SB')

	if ffs:
		logger.debug(f"Custom force fields provided: {ffs}")
		# ffs = ffs.split(",")
	else:
		logger.debug(f"No custom force fields provided. Using default: {ff_default}")
		ffs = ff_default

	# Prepare the complex (make sure the hydrogens naming matches amber)

	noH_complex = remove_hydrogens(pdbfile, logger_name, logfile)
	logger.debug(f"Hydrogens removed from {pdbpath}, saved as {noH_complex}")
	protonated_complex = add_hydrogens(noH_complex, logger_name, logfile, ffs)

	# Extract chain ids from the PDB file
	u = mda.Universe(pdbfile)
	chains = []
	for atom in u.atoms:
		if atom.segid and atom.segid not in chains:
			chains.append(atom.segid)
	logger.debug(f"Extracted chain IDs: {chains}")

	# Generate the receptor & ligand (also match the hydrogens to amber specifics)
	chain_files = gen_chains(noH_complex, logger_name, logfile, chains, ffs)
	solvate_gen_parm(protonated_complex, logger_name, logfile,chain_files, ffs)

	protonated_complex.replace('.pdb', '_solvated.prmtop')
	protonated_complex.replace('.pdb', '_solvated.pdb')

	end_time = time.time()
	elapsed_time = end_time - start_time
	logger.debug(f"Total execution time: {elapsed_time:.2f} seconds")

	# Check output files
	logger.debug("Checking generated files.")

	missing = [p for p in expected_files if not p.exists()]
	for file in expected_files:
		if file not in missing:
			logger.info(f"generated {file}")
		else:
			logger.error(f"{file} not generated")

	if not missing:
		checkpoint  = Path(f"check_01_amber_{int(time.time())}.done")
		checkpoint.parent.mkdir(exist_ok=True)
		with checkpoint.open('w') as f:
			msg = f"AmberTools preparation completed successfully for structure {structure}.\n"
			f.write(msg)
			logger.info(msg)
			os.chdir(cwd)
		return
	else:
		logger.error("AmberTools preparation failed. There are missing files.")
		sys.exit(1)
			
	# Return to initial working directory
	logger.debug(f"The current directory is {Path.cwd()}")
	os.chdir(cwd)
	logger.debug(f"Changed directory to {Path.cwd()}")


if __name__ == "__main__":
	main()


