#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import subprocess
import time
import sys
from pathlib import Path
import MDAnalysis as mda
from scripts.utils.helpers import setup_logger
from scripts.utils.forcefields import select_forcefields

# Number of tleap error lines repeated at ERROR level when tleap fails.
TLEAP_MAX_ERRORS = 20


def remove_hydrogens(filename: str, logger_name: str, logfile: Path):
	""" Remove hydrogens from the structure and prepare it with pdb4amber. """
	logger = setup_logger(logger_name, logfile)
	logger.debug(f"Removing hydrogens from {filename} using pdb4amber.")
	outfile = filename.replace('.pdb', '_noH.pdb')

	cmd = f"""
	pdb4amber -i {filename} -o {outfile} --nohyd --no-conect
    """
	
	process = subprocess.run(cmd, shell=True, text=True, capture_output=True)
	if process.returncode != 0:
		logger.error(f"There was an error removing hydrogens from {filename} to create {outfile}. Check the logs for the comple pdb4amber error message.")
		logger.debug(f"pdb4amber error message: \n {process.stderr}")
		return None  
	else:
		logger.debug(f"Hydrogens removed with pdb4amber. {outfile} was created.")
		logger.debug(f"pdb4amber log: \n {process.stdout}")
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
	# Errors are logged by _run_tleap; the caller checks that the output file exists.
	if _run_tleap(tleap_input, logger, f"adding hydrogens to {file}"):
		logger.debug(f"Hydrogens added with tleap. {outfile} was created.")

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
		chain_file = complex_name + f'_{chain}_noH.pdb'

		# Run pdb_selchain.py command for chain
		cmd_chain = ['pdb_selchain', f'-{chain}', complex_file]
		process_chain = subprocess.run(cmd_chain, capture_output=True, text=True)

		# Check for errors
		if process_chain.returncode != 0:
			logger.error(f"Error occurred while generating receptor file:\n {process_chain.stderr}")
			return None
		else:
			logger.debug(f"Chain {chain} extraction complete.")

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

	# True if tleap built and saved the systems without errors.
	return _run_tleap(cmd, logger, f"solvating {complex_file}")


def _run_tleap(cmd: str, logger, step: str) -> bool:
	"""
	Run tleap with the given input and log its output.

	The full output is logged at DEBUG (runs of more than 10 identical lines are
	collapsed). tleap's error lines (starting with Error, ERROR or FATAL) and its
	"Exiting LEaP: Errors = N" summary are logged at ERROR, the first TLEAP_MAX_ERRORS
	of them followed by a count of the rest. Returns False if tleap exits with an
	error code or reports errors, True otherwise.
	"""
	process = subprocess.run(['tleap', '-f', '-'], input=cmd, capture_output=True, text=True)
	lines = (process.stdout + process.stderr).splitlines()

	# Full output at DEBUG, with long runs of identical lines collapsed.
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
	logger.debug(f"tleap output ({step}):\n" + '\n'.join(filtered))

	# Error count from the summary line, e.g. "Exiting LEaP: Errors = 1; Warnings = 2; Notes = 1."
	summary = next((line.strip() for line in lines if line.strip().startswith("Exiting LEaP:")), "")
	match = re.search(r"Errors\s*=\s*(\d+)", summary)
	n_errors = int(match.group(1)) if match else 0
	if process.returncode == 0 and n_errors == 0:
		return True

	# tleap failed: its error lines at ERROR, the first ones and a count of the rest.
	errors = [line.strip() for line in lines if line.strip().startswith(("Error", "ERROR", "FATAL"))]
	logger.error(f"tleap failed while {step} (exit code {process.returncode}). {summary}")
	for line in errors[:TLEAP_MAX_ERRORS]:
		logger.error(f"tleap: {line}")
	if len(errors) > TLEAP_MAX_ERRORS:
		logger.error(f"tleap: ... and {len(errors) - TLEAP_MAX_ERRORS} more error lines (see the DEBUG output or leap.log).")
	return False



def _stop(logger, cwd: Path, message: str):
	"""Log an error, return to the starting directory and stop the pipeline."""
	logger.error(message)
	os.chdir(cwd)
	sys.exit(1)


########## MAIN ##########

def main(logfile: Path = None, pdbpath: Path = None,
	dna : bool = False, phos : bool = False, ffs: list = None,
	protein : bool = True, rna : bool = False):

	"""
    Prepare system using AmberTools. Add solvent and ions, generate prmtop
    and inpcrd files based on force fields.

    Arguments:
        logfile	    : Path to the log file.
        path-pdb    : Input pdb file path
        protein     : True if the system contains protein. Adds protein.ff19SB force field.
        phos        : True if the system contains phosphorylated residues. Adds phosaa19SB force field.
        dna         : True if the system contains DNA. Adds DNA.bsc1 force field.
        rna         : True if the system contains RNA. Adds RNA.OL3 force field.
        ffs         : List of force fields to use instead of the automatic choice
                        (names as in LEaP, e.g. protein.ff19SB). Avoid using clashing
                        force fields. water.opc is added if no water model is listed.
    The water model (water.opc) is always loaded, after the other force fields.
    """

	start_time = time.time()

	##### Initialize logger
	logger_name = 'prep'
	logger = setup_logger(logger_name, logfile)

	amberdir = Path(pdbpath).parent
	pdbfile = Path(pdbpath).name
	structure = Path(pdbpath).stem

	logger.debug(f"Extracted amberdir: {amberdir}, pdbfile: {pdbfile}, structure: {structure}")


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
		logger.debug(f"Checkpoint file(s) found: {', '.join(f.name for f in checkpoint_files)}")
		# Check if expected output files exist
		missing = [p for p in expected_files if not p.exists()]
		if not missing:
			logger.warning("The expected output files are already present. Skipping AmberTools preparation.")
			os.chdir(cwd)
			return

	# Check if input file is a PDB file
	if not pdbpath.suffix.lower() == '.pdb':
		logger.error(f"Input file must be a PDB file with .pdb extension. Check that {pdbpath} meets this criteria.")
		sys.exit(1)

	# Force fields for the molecule types present, water model last.
	logger.debug(f"Molecule types: protein={protein}, phos={phos}, dna={dna}, rna={rna}")
	if ffs:
		logger.info(f"Custom force fields provided: {ffs}")
	ffs = select_forcefields(protein=protein, phos=phos, dna=dna, rna=rna, ffs=ffs, logger=logger)
	logger.info(f"Force fields loaded in LEaP: {ffs}")

	# Prepare the complex (make sure the hydrogens naming matches amber)
	noH_complex = remove_hydrogens(pdbfile, logger_name, logfile)
	# pdb4amber failed: its error is already logged, the later steps need its output.
	if noH_complex is None or not Path(noH_complex).is_file():
		_stop(logger, cwd, f"pdb4amber could not process {pdbfile}. Check the structure for missing residues or atoms before running prep again.")
	protonated_complex = add_hydrogens(noH_complex, logger_name, logfile, ffs)
	# tleap can exit normally without writing the PDB, so the output file is checked.
	if not Path(protonated_complex).is_file():
		_stop(logger, cwd, f"tleap did not write {protonated_complex}. Check the log files for the comple tleap output.")

	# Extract chain ids from the PDB file
	u = mda.Universe(pdbfile)
	chains = []
	for atom in u.atoms:
		if atom.segid and atom.segid not in chains:
			chains.append(atom.segid)
	logger.debug(f"Extracted chain IDs: {chains}")

	# Generate the receptor & ligand (also match the hydrogens to amber specifics)
	chain_files = gen_chains(noH_complex, logger_name, logfile, chains, ffs)
	# Every chain needs its own protonated file for the LEaP step below.
	if chain_files is None:
		_stop(logger, cwd, "Chain extraction failed. See the error above.")
	missing_chains = [f for f in chain_files if not Path(f).is_file()]
	if missing_chains:
		_stop(logger, cwd, f"tleap did not write the chain files: {', '.join(missing_chains)}")
	# tleap's errors are already logged; the working directory is restored before stopping.
	if not solvate_gen_parm(protonated_complex, logger_name, logfile, chain_files, ffs):
		_stop(logger, cwd, f"tleap could not build the solvated system. See the errors above and leap.log in {amberdir}.")

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
			logger.debug(f"generated {file}")
		else:
			logger.error(f"{file} not generated")

	if not missing:
		checkpoint  = Path(f"check_01_amber_{int(time.time())}.done")
		checkpoint.parent.mkdir(exist_ok=True)
		with checkpoint.open('w') as f:
			msg = f"AmberTools preparation completed successfully for structure {structure}."
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


