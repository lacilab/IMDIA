#!/bin/bash

source "../utils/colors.sh"

# Print usage/help
print_usage() {
    cat <<'USAGE'
Usage: wrap.sh gromacs_directory
Refine a wrapped and fitted trajectory from the MD simulation.

Positional arguments:
  gromacs_directory   Directory containing GROMACS input files (md.tpr, md.xtc, index.ndx etc.)

Options:
  --with-dna          Complex includes DNA
  --with-phosphate    Complex includes phosphorylated residues
  -h, --help          Show this help message and exit
USAGE
}

# Check args
if [ "$#" -ne 1 ]; then
    print_usage
    exit 1
fi

cd $1

# LOG
LOGFILE="wrap_$rep{n}_$(date +%Y%m%d_%H%M%S).log"
# Redirect all stdout and stderr to tee
# This will print to screen AND save to the logfile, removing color tags
exec > >(tee >(sed -r "s/\x1B\[[0-9;]*[mK]//g" >> "$LOGFILE")) 2>&1
echo "Logging to: $LOGFILE"

# Check for checkpoint files
if [ -f "check_03_wrap_*.done" ]; then
    warn "Checkpoint file found. Skipping wrapping."
    echo "To rerun this step, please remove the checkpoint file(s): check_03_wrap_*.done"
    exit 0
fi

# Default options
with_dna=false
with_phosphate=false

# Parse flags before positional arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-dna)
      with_dna=true
      echo "Including --with-dna flag. System should include DNA."
      shift
      ;;
    --with-phosphate)
      with_phosphate=true
      echo "Including --with-phosphate flag. System should include phosphorylated residues."
      shift
      ;;
      ;;
    -h|--help)
      print_usage
      exit 0
      ;;
    -*)
      echo "Unknown option: $1"
      print_usage
      exit 1
      ;;
    *)
      break
      ;;
  esac
done

if [[ $with_dna == true ]] && [[ $with_phosphate == true ]]; then
    protein_selection="Protein_DNA_PTR"
else
    protein_selection="Protein"
fi

# center molecule and put all atoms at the closest distance from center of the box
# TODO check manually with active and inactive structures, added -center from the previous version
gmx trjconv -s md.tpr -f md.xtc -o md_noPBC_mol_com.xtc -pbc mol -ur compact -center -quiet -n index.ndx << EOL
0
EOL

# Cluster to avoid monomers jumping 
gmx trjconv -s md.tpr -f md_noPBC_mol_com.xtc -o md_noPBC_mol_com_cluster.xtc -pbc cluster -quiet -n index.ndx << EOL
$protein_selection
0
EOL

# Fit the system to reference structure (using Backbone)
gmx trjconv -s md.tpr -f md_noPBC_mol_com_cluster.xtc -o md_noPBC_mol_com_cluster_fit.xtc -fit rot+trans -quiet -n index.ndx << EOL
Backbone
0
EOL

# Clean up intermediate files
rm md_noPBC_mol_com.xtc md_noPBC_mol_com_cluster.xtc

# Output .pdb files for original trajectory reference (t=0 of raw simulation)
gmx trjconv -s md.tpr -f md.xtc -o md_sys.pdb -pbc mol -ur compact -dump 0 -quiet -n index.ndx << EOL
0
EOL

gmx trjconv -s md.tpr -f md.xtc -o md_prot.pdb -pbc mol -ur compact -dump 0 -quiet -n index.ndx << EOL
$protein_selection
EOL

# Generate checkpoint file
touch "check_03_wrap_$(date +%Y%m%d_%H%M%S).done"

