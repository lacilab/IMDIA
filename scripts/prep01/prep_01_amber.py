#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import subprocess
import time
import sys
import math
from pathlib import Path
import MDAnalysis as mda
from scripts.utils.helpers import setup_logger
from scripts.prep01.prep_01b_protonation import (
    assign_protonation, canonical, is_hydrogen, residue_key, residues, xyz,
    terminal_residues, PROTEIN_NAMES, PHOSPHO_RESIDUES, NUCLEIC_NAMES,
    WATER_NAMES, nucleic_base,
)
from scripts.utils.forcefields import select_forcefields
from scripts.utils import ligands

# Number of tleap error lines repeated at ERROR level when tleap fails.
TLEAP_MAX_ERRORS = 20

# Helper to match residue IDs by order original and prepared

def _residue_pairs(reference, prepared):
    """Match ordered residues despite changes in numbering or chain IDs."""
    left, right = residues(reference), residues(prepared)

    if len(left) != len(right):
        raise ValueError(
            'Preparation changed residue count; explicit mapping is required'
        )

    pairs = []
    for (key, name), (other, converted) in zip(left.items(), right.items()):
        matching = canonical(name) == canonical(converted)

        if name in WATER_NAMES and converted in WATER_NAMES:
            matching = True

        if name in NUCLEIC_NAMES and converted in NUCLEIC_NAMES:
            matching = nucleic_base(name) == nucleic_base(converted)

        if not matching:
            raise ValueError(
                f'Preparation changed residue identity/order: '
                f'{key} {name} -> {converted}'
            )

        pairs.append((key, other))

    return pairs

# Helper to remove newly generated TER records
def _restore_peptide_continuity(original, prepared, logger):
    # Remove newly added TER records only where the original segment was
    # continuous and the prepared backbone C–N distance supports a peptide bond.
    source = Path(original).read_text().splitlines(keepends=True)
    target = Path(prepared).read_text().splitlines(keepends=True)
    mapping = {new: old for old, new in _residue_pairs(source, target)}

    segments, segment, previous_chain = {}, 0, None
    for line in source:
        if line.startswith('TER'):
            segment += 1
            previous_chain = None
        elif line.startswith(('ATOM  ', 'HETATM')):
            key = residue_key(line)
            if previous_chain is not None and key[0] != previous_chain:
                segment += 1
            segments[key], previous_chain = segment, key[0]

    atoms, following, next_key = {}, {}, None
    for i in range(len(target) - 1, -1, -1):
        line = target[i]
        following[i] = next_key
        if line.startswith(('ATOM  ', 'HETATM')):
            next_key = residue_key(line)
            atoms[next_key, line[12:16].strip()] = line

    names, previous, removed, output = residues(target), None, 0, []
    for i, line in enumerate(target):
        if line.startswith(('ATOM  ', 'HETATM')):
            previous = residue_key(line)

        if (
            line.startswith('TER')
            and previous is not None
            and following[i] is not None
        ):
            next_key = following[i]
            protein = all(
                canonical(names[key]) in PROTEIN_NAMES
                or names[key] in PHOSPHO_RESIDUES
                for key in (previous, next_key)
            )

            if (
                protein
                and previous[0] == next_key[0]
                and segments[mapping[previous]] == segments[mapping[next_key]]
            ):
                carbon = atoms.get((previous, 'C'))
                nitrogen = atoms.get((next_key, 'N'))

                if (
                    carbon and nitrogen
                    and 1.0 <= math.dist(xyz(carbon), xyz(nitrogen)) <= 1.8
                ):
                    removed += 1
                    continue

        output.append(line)

    if removed:
        path = Path(prepared)
        path.with_suffix('.before_TER_cleanup.pdb').write_text(''.join(target))
        path.write_text(''.join(output))
        logger.info(
            f'Removed {removed} added TER records across continuous peptide '
            'bonds; original breaks retained.'
        )

def _complete_terminal_oxygens(prepared, built, logger):
	# Copy verified missing terminal OXT atoms from LEaP, preserving existing atoms.
	path = Path(prepared)
	lines = path.read_text().splitlines(keepends=True)
	donor = Path(built).read_text().splitlines(keepends=True)
	pairs = _residue_pairs(lines, donor)

	atoms, last, built_atoms = {}, {}, {}
	for i, line in enumerate(lines):
		if line.startswith(('ATOM  ', 'HETATM')):
			key = residue_key(line)
			atom_key = key, line[12:16].strip()
			if atom_key in atoms:
				raise ValueError(f'Duplicate prepared atom: {atom_key}')
			atoms[atom_key], last[key] = line, i

	for line in donor:
		if line.startswith(('ATOM  ', 'HETATM')) and not is_hydrogen(line):
			atom_key = residue_key(line), line[12:16].strip()
			if atom_key in built_atoms:
				raise ValueError(f'Duplicate LEaP atom: {atom_key}')
			built_atoms[atom_key] = line

	endpoints, names = terminal_residues(lines), residues(lines)
	nitrogens = [
		(key, xyz(line))
		for (key, atom), line in atoms.items()
		if atom == 'N'
	]
	serial = max((int(line[6:11]) for line in atoms.values()), default=0)
	additions = {}

	for key, other in pairs:
		oxygen = built_atoms.get((other, 'OXT'))
		if oxygen is None or (key, 'OXT') in atoms:
			continue

		if (
			canonical(names[key]) not in PROTEIN_NAMES
			and names[key] not in PHOSPHO_RESIDUES
		):
			continue

		carbon = atoms.get((key, 'C'))
		if key not in endpoints or carbon is None or (key, 'O') not in atoms:
			raise ValueError(f'Unexpected LEaP terminal oxygen at {key}')

		if any(
			other_key != key and math.dist(xyz(carbon), point) <= 1.8
			for other_key, point in nitrogens
		):
			raise ValueError(
				f'LEaP added OXT to a peptide-bonded carbonyl at {key}; '
				'review chain breaks'
			)

		built_carbon = built_atoms.get((other, 'C'))
		if (
			built_carbon is None
			or math.dist(xyz(carbon), xyz(built_carbon)) > 0.002
			or not 1.0 <= math.dist(xyz(carbon), xyz(oxygen)) <= 1.6
		):
			raise ValueError(
				f'Invalid LEaP terminal oxygen coordinates at {key}'
			)

		serial += 1
		if serial > 99999:
			raise ValueError(
				'Terminal completion exceeds the PDB atom-serial limit'
			)

		additions[last[key]] = (
			oxygen[:6] + f'{serial:5d}' + oxygen[11:17]
			+ carbon[17:27] + oxygen[27:]
		)

	if additions:
		path.with_suffix('.before_terminal_completion.pdb').write_text(
			''.join(lines)
		)
		path.write_text(
			''.join(line + additions.get(i, '') for i, line in enumerate(lines))
		)
		logger.info(
			f'Completed {len(additions)} missing terminal oxygen atoms '
			'using the loaded LEaP templates.'
		)

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
	    # Repair added chain breaks before passing the structure to LEaP.
	try:
		_restore_peptide_continuity(filename, outfile, logger)
	except (ValueError, OSError) as exc:
		logger.error(f'Cannot establish peptide continuity: {exc}')
		return None
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
		check complex
charge complex
		savePdb complex {outfile}
		quit
	"""
	logger.debug(f"Running tleap with input:\n{tleap_input}")
	# Errors are logged by _run_tleap; the caller checks that the output file exists.
	ok, lines = _run_tleap(tleap_input, logger, f"adding hydrogens to {file}")
	if ok:
		logger.debug(f"Hydrogens added with tleap. {outfile} was created.")
	# Residues LEaP has no template for (they load, but get no atom types).
	return (outfile if ok else None), _unknown_residues(lines)


def _unknown_residues(lines: list[str]) -> list[tuple[str, int]]:
	"""
	Residues tleap reported as having no template while loading a PDB file.

	Returns (residue name, LEaP sequence number) pairs, unique and in order, from
	"Creating new UNIT for residue: <NAME> sequence: <N>" lines, or from
	"Unknown residue: <NAME> number: <N>" lines if the first kind is absent.
	"""
	found = []
	for pattern in (r"Creating new UNIT for residue:\s*(\S+)\s+sequence:\s*(\d+)",
					r"Unknown residue:\s*(\S+)\s+number:\s*(\d+)"):
		for line in lines:
			match = re.search(pattern, line)
			if match and (match.group(1), int(match.group(2))) not in found:
				found.append((match.group(1), int(match.group(2))))
		if found:
			break
	return found


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

		chain_file, _ = add_hydrogens(chain_file, logger_name, logfile, ffs)
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
	ok, _ = _run_tleap(cmd, logger, f"solvating {complex_file}")
	return ok


def _run_tleap(cmd: str, logger, step: str) -> tuple[bool, list[str]]:
	"""
	Run tleap with the given input and log its output.

	The full output is logged at DEBUG (runs of more than 10 identical lines are
	collapsed). tleap's error lines (starting with Error, ERROR or FATAL) and its
	"Exiting LEaP: Errors = N" summary are logged at ERROR, the first TLEAP_MAX_ERRORS
	of them followed by a count of the rest. Returns (ok, output lines); ok is False
	if tleap exits with an error code or reports errors.
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
		return True, lines

	# tleap failed: its error lines at ERROR, the first ones and a count of the rest.
	errors = [line.strip() for line in lines if line.strip().startswith(("Error", "ERROR", "FATAL"))]
	logger.error(f"tleap failed while {step} (exit code {process.returncode}). {summary}")
	for line in errors[:TLEAP_MAX_ERRORS]:
		logger.error(f"tleap: {line}")
	if len(errors) > TLEAP_MAX_ERRORS:
		logger.error(f"tleap: ... and {len(errors) - TLEAP_MAX_ERRORS} more error lines (see the DEBUG output or leap.log).")
	return False, lines



def _stop(logger, cwd: Path, message: str):
	"""Log an error, return to the starting directory and stop the pipeline."""
	logger.error(message)
	os.chdir(cwd)
	sys.exit(1)


########## MAIN ##########

def main(logfile: Path = None, pdbpath: Path = None,
	dna : bool = False, phos : bool = False, ffs: list = None,
	protein : bool = True, rna : bool = False, protonation: dict = None):

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
        protonation : Protonation configuration (default: None)
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

	# Remove hydrogens and prepare the PDB.
	noH_complex = remove_hydrogens(pdbfile, logger_name, logfile)
	# pdb4amber failed: its error is already logged, the later steps need its output.
	if noH_complex is None or not Path(noH_complex).is_file():
		_stop(logger, cwd, f"pdb4amber could not process {pdbfile}. Check the structure for missing residues or atoms before running prep again.")
	protonated_complex, unknown = add_hydrogens(noH_complex, logger_name, logfile, ffs)

	# Residues without a template in the loaded force fields are removed, then hydrogens are added again.
	removed = []
	if unknown:
		names = sorted({name for name, _ in unknown})
		# Residues inside a protein or nucleic-acid chain are not removed (it would break the chain).
		in_chain = ligands.polymer_residues(noH_complex, names)
		if in_chain:
			_stop(logger, cwd, f"Residues without a template in the loaded force fields are part of a "
				  f"polymer chain: {', '.join(in_chain)}. Add a force field or parameters for them.")
		# The unfiltered file is kept; the filtered one takes its name, so later file names are unchanged.
		unfiltered = noH_complex.replace('_noH.pdb', '_noH_unfiltered.pdb')
		os.replace(noH_complex, unfiltered)
		removed = ligands.remove_residues(unfiltered, names, noH_complex)
		logger.warning(f"VERY IMPORTANT. IMDIA removed these residues from the structure because they do not have a template in the loaded force fields: {', '.join(removed)}. "
					   f"\nThe structure before removal is {unfiltered}."
						"\nKeeping this residues with the current settings would have caused the pipeline to fail, so they were removed by default."
						"\nIf you want to run the simulation with these residues, include a forcefield that supports them in the parameters."
						"\nAlternatively, you can try preparing the system with a suitable package or web based tool, and using the prepared files to run the MD simulation in IMDIA."
				   "\nAdd parameters for them if they are needed in the simulation. ")
		protonated_complex, unknown = add_hydrogens(noH_complex, logger_name, logfile, ffs)
		if unknown:
			_stop(logger, cwd, "tleap still reports residues without a template after removal: "
				  + ", ".join(f"{name} {number}" for name, number in unknown))
	# Keep the reference consistent if unsupported ligands were removed.
	reference_pdb = pdbfile
	if removed:
		reference_pdb = f"{structure}_filtered_original.pdb"
		ligands.remove_residues(pdbfile, names, reference_pdb)

	# Complete terminal heavy atoms before PROPKA and PDB2PQR read the structure.
	if protonated_complex is None or not Path(protonated_complex).is_file():
		_stop(
			logger, cwd,
			'LEaP could not complete the initial structure; inspect leap.log.',
		)

	try:
		_complete_terminal_oxygens(noH_complex, protonated_complex, logger)
	except (ValueError, OSError) as exc:
		_stop(logger, cwd, f'Terminal oxygen completion failed: {exc}')
	# Assign protonation after ligand filtering.
	assigned_complex = assign_protonation(
		pdbfile=noH_complex,
		original_pdb=reference_pdb,
		config=protonation,
		logger=logger,
	)

	# Build hydrogens using the assigned states.
	protonated_complex, unknown = add_hydrogens(
		assigned_complex, logger_name, logfile, ffs
	)
	if unknown:
		_stop(
			logger, cwd,
			"Assigned residues have no force-field template: "
			+ ", ".join(
				f"{name} {number}" for name, number in unknown
			)
		)
	# tleap can exit normally without writing the PDB, so the output file is checked.
	if protonated_complex is None or not Path(protonated_complex).is_file():
		_stop(logger, cwd, f"tleap did not write {protonated_complex}. Check the log files for the comple tleap output.")

	# Extract chain IDs from the assigned structure.
	u = mda.Universe(assigned_complex)
	chains = []
	for atom in u.atoms:
		if atom.segid and atom.segid not in chains:
			chains.append(atom.segid)
	# A chain whose residues were all removed gets no chain file.
	if removed:
		present = {line[21] for line in Path(noH_complex).read_text().splitlines()
				   if line.startswith(("ATOM", "HETATM"))}
		chains = [c for c in chains if len(c) != 1 or c in present]
	logger.debug(f"Extracted chain IDs: {chains}")

	# Generate the receptor & ligand (also match the hydrogens to amber specifics)
	chain_files = gen_chains(assigned_complex, logger_name, logfile, chains, ffs)
	# Every chain needs its own protonated file for the LEaP step below.
	if chain_files is None:
		_stop(logger, cwd, "Chain extraction failed. See the error above.")
	missing_chains = [str(f) for f in chain_files if f is None or not Path(f).is_file()]
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


