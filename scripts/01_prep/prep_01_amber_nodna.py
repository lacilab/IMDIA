import os
import subprocess
import time
import sys

def remove_hydrogens(file):
    cmd = f'''
pdb4amber -i {file}.pdb -o {file}_noH.pdb --nohyd --no-conect
pdb4amber -i {file}_noH.pdb -o {file}_reduced.pdb --reduce --no-conect
    '''

    process = subprocess.run(cmd, shell=True, text=True, capture_output=True)

    if process.returncode != 0:
        print("Error:", process.stderr)
        return None  
    else:
        print(f"Hydrogen processing complete. {process.stdout}")

    return f"{file}_reduced.pdb"


def add_hydrogens(file):
	"""Add hydrogen atoms to a PDB file using tleap."""
	protonated_file = file.replace('_reduced.pdb', '_protonated.pdb')
	tleap_input = f"""
source leaprc.protein.ff19SB
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
		mono_a='A', mono_b='B'
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


def solvate_gen_parm(complex_file, mono_a_file, mono_b_file):
	"""Protonate and solvate the complex and chains using LEaP."""
	complex_name = os.path.splitext(complex_file)[0]
	complex_solvated_file = complex_name + '_solvated.pdb'

	mono_a_filebase = mono_a_file[:-4]
	mono_b_filebase = mono_b_file[:-4]

	cmd = f'''
source leaprc.protein.ff19SB
source leaprc.water.opc
loadamberparams frcmod.ionslm_1264_opc

complex = loadPdb {complex_file}
mono_a = loadPdb {mono_a_file}
mono_b = loadPdb {mono_b_file}

set default PBRadii mbondi2

saveAmberParm complex {complex_name}.prmtop {complex_name}.inpcrd
saveAmberParm mono_a {mono_a_filebase}.prmtop {mono_a_filebase}.inpcrd
saveAmberParm mono_b {mono_b_filebase}.prmtop {mono_b_filebase}.inpcrd

solvateOct complex OPCBOX 10.0
# Add ions
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

# struct_dir = sys.argv[1]
# complex_file = f'{struct_dir}'
# path_complex = os.path.join(struct_dir, complex_file)

path_complex = sys.argv[1]

# Prepare the complex (make sure the hydrogens naming matches amber)
noH_complex = remove_hydrogens(path_complex)
protonated_complex = add_hydrogens(noH_complex)

# Generate the receptor & ligand (also match the hydrogens to amber specifics)
mono_a_file, mono_b_file = gen_chains(noH_complex, 'A', 'B')
solvate_gen_parm(protonated_complex, mono_a_file, mono_b_file)
solvated_complex_prmtop = protonated_complex.replace('.pdb', '_solvated.prmtop')
solvated_complex_pdb = protonated_complex.replace('.pdb', '_solvated.pdb')

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Total execution time: {elapsed_time:.2f} seconds")
