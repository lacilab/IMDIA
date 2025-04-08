#!/bin/bash
source /opt/gromacs/gromacs_2024/bin/GMXRC

# Check args
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 gromacs_directory"
    exit 1
fi

cd $1

# center all molecules in box and put all atoms at the closest distance from center of the box
gmx trjconv -s md.tpr -f md.xtc -o md_noPBC_mol_com.xtc -pbc mol -ur compact -quiet -n index.ndx << EOL
0
EOL

# Cluster to avoid monomers jumping 
gmx trjconv -s md.tpr -f md_noPBC_mol_com.xtc -o md_noPBC_mol_com_cluster.xtc -pbc cluster -quiet -n index.ndx << EOL
Protein_DNA_PTR
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
Protein_DNA_PTR
EOL

