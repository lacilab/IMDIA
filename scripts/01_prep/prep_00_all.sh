#!/bin/bash

source activate base
conda activate thesis_md
source /opt/gromacs/gromacs_2024/bin/GMXRC

# Check if datadir argument is provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 datadir filename"
    exit 1
fi

oridir=$(pwd)
datadir=$1
file=$2

# Remove lines starting with CONECT 
sed '/^CONECT/d' "$datadir/${file}.pdb" > "$datadir/${file}_noconect.pdb"

mkdir -p "$datadir/AmberTools"
mv "$datadir/${file}_noconect.pdb" "$datadir/AmberTools/"

########## PREPARATION ##########

# AmberTools preparation
cp prep_01_amber.py "$datadir/AmberTools/"
cd "$datadir/AmberTools/"
python prep_01_amber.py "${file}_noconect"
rm prep_01_amber.py
cd $oridir

# Convert to Gromacs format with parmed
mkdir "$datadir/gromacs_common"
gmxdir="$datadir/gromacs_common"
python prep_02_amb2gro.py -p "${datadir}/AmberTools/${file}_noconect_protonated_solvated.prmtop" \
    -c "$datadir/AmberTools/${file}_noconect_protonated_solvated.inpcrd" \
    -t "$gmxdir/topol.top" \
    -g "$gmxdir/gromacs.gro" \
#     -b "$gmxdir/gromacs.pdb"
# 
# Extract itps
mkdir "$gmxdir/toppar"
bash prep_03_extractitp.sh "$gmxdir"

# Remove the 2 ions added with AmberTools
bash prep_04_removeamberions.sh "$gmxdir"

# update topology file
sed -i "/Na/d" "$gmxdir/topol.top"
sed -i "/Cl/d" "$gmxdir/topol.top"

# Change ion names to gromacs standards
sed -i "s/Na+/Na/g" $gmxdir/toppar/molecule_5.itp $gmxdir/toppar/forcefield.itp 
sed -i "s/Cl-/Cl/g" $gmxdir/toppar/molecule_6.itp $gmxdir/toppar/forcefield.itp 

cp ../02_md/ions.mdp $gmxdir
 
# Add ions with GROMACS to neutralize and reach target concentration
gmx grompp -f "$gmxdir/ions.mdp" -c "$gmxdir/updated_gromacs.gro" -p "$gmxdir/topol.top" -o "$gmxdir/ions.tpr" -quiet -maxwarn 1
# Change the name of the solvent from WAT to Water in .top so it can be selected as a group.
sed -i "s/WAT/Water/g" $gmxdir/toppar/molecule_7.itp $gmxdir/toppar/forcefield.itp $gmxdir/topol.top
# There will be a warning because the system is not neutralized.
gmx genion -s "$gmxdir/ions.tpr" -o "$gmxdir/ions.gro" -p "$gmxdir/topol.top" -pname Na+ -nname Cl- -neutral -conc 0.15 -quiet << EOL
Water
quit
EOL

# Add DNA H-bonds constrains
cp prep_05_restraindnahbonds.sh $gmxdir
cd $gmxdir
bash prep_05_restraindnahbonds.sh
rm prep_05_restraindnahbonds.sh

# Change ion names to gromacs standards
sed -i "s/Cl-/Cl/g" topol.top
sed -i "s/Na+/Na/g" topol.top
sed -i "s/Water/WAT/g" toppar/molecule_7.itp toppar/forcefield.itp topol.top