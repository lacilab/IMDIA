import os
import subprocess
import time
import sys
import argparse

def remove_hydrogens(dir, filename):
    cmd = f'''
pdb4amber -i {dir}/{filename}.pdb -o {dir}/{filename}_noH.pdb --nohyd --no-conect
pdb4amber -i {dir}/{filename}_noH.pdb -o {dir}/{filename}_reduced.pdb --reduce --no-conect
    '''

    process = subprocess.run(cmd, shell=True, text=True, capture_output=True)

    if process.returncode != 0:
        print("Error:", process.stderr)
        return None  
    else:
        print(f"Hydrogen processing complete. {process.stdout}")

    return f"{dir}/{filename}_reduced.pdb"


def add_hydrogens(file, ffs=['protein.ff19SB', 'water.opc']):
	"""Add hydrogen atoms to a PDB file using tleap."""
	protonated_file = file.replace('_reduced.pdb', '_protonated.pdb')

	# Build tleap input dynamically
	ffs_import = "\n".join(f"source leaprc.{ff}" for ff in ffs)

	tleap_input = f"""
		{ffs_import}
		complex = loadPdb {file}
		savePdb complex {protonated_file}
		quit
	"""
	print(f"Running tleap with input:\n{tleap_input}")  # Debugging: Show tleap input
	process = subprocess.run(['tleap', '-f', '-'], input=tleap_input, text=True, capture_output=True)
	
	if process.returncode != 0:
		print(f"tleap error: {process.stderr}")
	else:
		print(f"tleap output: {process.stdout}")

	return protonated_file


def gen_chains(complex_file, *chains):
	"""Generate separate PDB files.
		*chains takes a varible number of arguments for the chains to be separated, for example:
		mono_a='A', mono_b='B', dna_a='C', dna_b='D'
	"""
	complex_name = os.path.basename(complex_file)

	files = []

	for chain in chains:
		chain_file = complex_file.replace(complex_name, f'{chain}_noH.pdb')

		# Run pdb_selchain.py command for chain
		cmd_chain = ['pdb_selchain', f'-{chain}', complex_file]
		process_chain = subprocess.Popen(cmd_chain, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		stdout_chain, stderr_chain = process_chain.communicate()

		# Check for errors
		if process_chain.returncode != 0:
			print(f"Error occurred while generating receptor file: {stderr_chain}")
			return None

		# Write receptor data to file
		with open(chain_file, 'w') as f:
			f.write(stdout_chain)

		chain_file = add_hydrogens(chain_file)
		files.append(chain_file)

	return tuple(files)


def solvate_gen_parm(complex_file, ffs=['protein.ff19SB', 'water.opc'], *chains):
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
		cmd += f'{chain} = loadPdb {chain}\n'

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




########## RUN ##########

start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('--path-pdb', default='.', help='Input pdb file path')
parser.add_argument('--with-dna', action='store_true', help='Use if the system contains DNA. Adds DNA.bsc1 force field.')
parser.add_argument('--with-phosphate', action='store_true', help='Use if the system contains phosphorylated residues. Adds phosaa19SB force field.')
parser.add_argument('--ffs', nargs='*', help='List of force fields to use from the list available in LEaP. Avoid using clashing force fields. Default: protein.ff19SB, water.opc')

# print help if no arguments are provided
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(0)

args = parser.parse_args()

amberdir, pdb = os.path.split(args.path_pdb)
gmxdir = os.path.dirname(amberdir)
filename = os.path.splitext(pdb)[0]

ff_default = ['water.opc']

# Check if input file is a PDB file
if not pdb.lower().endswith('.pdb'):
	print("Error: Input file must be a PDB file with .pdb extension.")
	sys.exit(1)

if args.with_dna:
	print("Flag --with-dna used. This is a DNA-protein complex.")
	ff_default.append('DNA.bsc1')

if args.with_phosphate:
	print("Flag --with-phosphate used. This system contains phosphorylated residues.")
	ff_default.append('phosaa19SB')
else:
	ff_default.append('protein.ff19SB')

if args.ffs:
	print(f"Custom force fields provided: {args.ffs}")
	ffs = args.ffs
else:
	print(f"No custom force fields provided. Using default: {ff_default}")
	ffs = ff_default

# Prepare the complex (make sure the hydrogens naming matches amber)
noH_complex = remove_hydrogens(amberdir, filename)
protonated_complex = add_hydrogens(noH_complex, ffs)

# Generate the receptor & ligand (also match the hydrogens to amber specifics)
if args.incl_dna:
	mono_a_file, mono_b_file, dna_a_file, dna_b_file = gen_chains(noH_complex, 'A', 'B', 'C', 'D')
	solvate_gen_parm(protonated_complex, ffs, mono_a_file, mono_b_file, dna_a_file, dna_b_file)
else:
	mono_a_file, mono_b_file = gen_chains(noH_complex, 'A', 'B')	
	solvate_gen_parm(protonated_complex, ffs, mono_a_file, mono_b_file)

solvated_complex_prmtop = protonated_complex.replace('.pdb', '_solvated.prmtop')
solvated_complex_pdb = protonated_complex.replace('.pdb', '_solvated.pdb')

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Total execution time: {elapsed_time:.2f} seconds")

# Check output files
expected_files = [
	f"{amberdir}/{filename}_noconect_protonated_solvated.inpcrd",
	f"{amberdir}/{filename}_noconect_protonated_solvated.prmtop",
	f"{gmxdir}/topol.top",
	f"{gmxdir}/gromacs.gro"
]

print("\nChecking generated files.")

missing = [p for p in expected_files if not os.path.isfile(p)]
for file in expected_files:
	if file not in missing:
		print(f"OK: generated {file}")
	else:
		print(f"ERROR: {file} not generated")

if not missing:
	print("OK: AmberTools preparation completed successfully.")
	sys.exit(0)
else:
	print("ERROR: AmberTools preparation failed. There are missing files.")
	sys.exit(1)
