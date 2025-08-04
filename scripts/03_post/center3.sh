#!/bin/bash
source /opt/gromacs/gromacs_2024/bin/GMXRC

# Check args
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 gromacs_directory"
    exit 1
fi

cd $1


#STEP 1: make molecules whole
gmx trjconv -s md.tpr -f md.xtc -o step1_whole.xtc -pbc whole -quiet -n index.ndx << EOL
0
EOL

# STEP 2: Clustering
gmx trjconv -s md.tpr -f step1_whole.xtc -o step2_cluster.xtc  -pbc cluster -quiet -n index.ndx << EOL
Protein_DNA_PTR
0
EOL

#STEP 3: No jumps
gmx trjconv -s md.tpr -f  step2_cluster.xtc -o step3_nojump.xtc -pbc nojump -quiet -n index.ndx << EOL
0
EOL


#STEP 3: center DNA in box
gmx trjconv -s md.tpr -f md.xtc -o step1_centerdna.xtc -center -pbc mol -quiet -n index.ndx << EOL
DNA
0
EOL

#STEP 3: center complex in box
gmx trjconv -s md.tpr -f step1_centerdna.xtc -o step2_centercomplex.xtc -center -quiet -n index.ndx << EOL
Protein_DNA_PTR
0
EOL


#STEP 4: fit the system to reference in structure file
gmx trjconv -s md.tpr -f step2_centercomplex.xtc  -o step7_reference.xtc -fit rot+trans -quiet -n index.ndx << EOL
Backbone
0
EOL


#STEP 7: output a .pdb file that contains the exact particles as the trajectory at t=0, so other programs know how to interpret the numbers in the trajectory
gmx trjconv -s md.tpr -f md.xtc -o md_sys_fixed.pdb -pbc mol -ur compact -dump 0 -quiet -n index.ndx << EOL
0
EOL



