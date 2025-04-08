#!/bin/bash
source /opt/gromacs/gromacs_2024/bin/GMXRC

# Check if the top file is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 gromacs_directory replicate_n"
    exit 1
fi

scriptdir=$(pwd)
gmxdir=$1
n=$2

cd "$gmxdir"
cd ..
mkdir "rep$n"
mkdir "rep$n/gromacs"
cd "rep$n/gromacs"
outdir=$(pwd)

cd "$scriptdir"
cp *.mdp "$outdir"
cd "$outdir"

# Energy minimization
gmx grompp -f minim.mdp -c "$gmxdir/ions.gro" -r "$gmxdir/ions.gro" -p "$gmxdir/topol.top" -o em.tpr -quiet -maxwarn 1
gmx mdrun -v -deffnm em -quiet -nb gpu

# Create groups
gmx make_ndx -f em.gro -o index.ndx -quiet << EOL
r DA | r DC | r DT | r DG | r DC3 | r DC5 | r DG3 | r DG5
name 19 DNA
1 | 19 | 13
name 20 Protein_DNA_PTR
14 | 15 | 16
name 21 Water_and_ions 
quit
EOL

# NVT
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p "$gmxdir/topol.top" -o nvt.tpr -n index.ndx -quiet -maxwarn 1
gmx mdrun -v -deffnm nvt -quiet -nb gpu

# NPT
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p "$gmxdir/topol.top" -o npt.tpr -n index.ndx -quiet -maxwarn 1 
gmx mdrun -v -deffnm npt -quiet -nb gpu

# MD
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p "$gmxdir/topol.top" -o md.tpr -n index.ndx -quiet -maxwarn 1
gmx mdrun -ntmpi 16 -ntomp 8 -nb gpu -pme gpu -npme 1 -gputasks 0000000000111111 -v -deffnm md -quiet > output1.log 2>&1 &
# gmx mdrun -ntmpi 16 -ntomp 8 -cpi md.cpt -append -nb gpu -pme gpu -npme 1 -gputasks 0000000000111111 -v -deffnm md -quiet > output2.log 2>&1 &
