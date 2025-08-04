#!/bin/bash
source /opt/gromacs/gromacs_2024/bin/GMXRC

# Check args
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 gromacs_directory"
    exit 1
fi

cd $1



#STEP 1: center DNA in box
gmx trjconv -s md.tpr -f md.xtc -o step1_centerdna.xtc -center -pbc mol -quiet -n index.ndx << EOL
DNA
0
EOL

#STEP 2: center complex in box
gmx trjconv -s md.tpr -f step1_centerdna.xtc -o step2_centercomplex.xtc -center -pbc mol -ur compact -quiet -n index.ndx << EOL
Protein_DNA_PTR
0
EOL


#STEP 3: fit the system to reference in structure file
gmx trjconv -s md.tpr -f step2_centercomplex.xtc  -o step3_reference.xtc -fit rot+trans -quiet -n index.ndx << EOL
Backbone
0
EOL


#STEP 4: output a .pdb file that contains the exact particles as the trajectory at t=0, so other programs know how to interpret the numbers in the trajectory
gmx trjconv -s md.tpr -f md.xtc -o md_sys_fixed_rep.pdb -pbc mol -ur compact -dump 0 -quiet -n index.ndx << EOL
0
EOL



