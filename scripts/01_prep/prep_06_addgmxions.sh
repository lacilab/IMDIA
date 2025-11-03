#!/bin/bash

# Print usage/help
print_usage() {
    cat <<'USAGE'
Usage: prep_06_addgmxions.sh topology
Add GROMACS ions to the system to neutralize it and reach target concentration.

Positional arguments:
  gmxdir    Path to the GROMACS directory (containing gromacs.gro and topol.top)
  ions.mdp  Path to the GROMACS parameter file to add ions.

Options:
  -h, --help          Show this help message and exit
USAGE
}

# Check for help argument
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    print_usage
    exit 0
fi

# Check if the top file is provided
if [ "$#" -ne 1 ]; then
    echo "ERROR: gmxdir argument is required."
    print_usage
    exit 1
fi

gmxdir=$1
# default ions.mdp path
ions_mdp="../02_md/ions.mdp"
# assign second argument if provided
if [ "$#" -eq 2 ]; then
    ions_mdp="$2"
    # check if the provided ions.mdp exists and is an .mdp file
    if [ ! -f "$ions_mdp" ] || [[ "$ions_mdp" != *.mdp ]]; then
        echo "ERROR: Provided parameter file '$ions_mdp' does not exist or is not an .mdp file."
        exit 1
fi

cp "$ions_mdp" "$gmxdir"
# Check the ions.mdp is in the target folder
if [ -f "$gmxdir/ions.mdp" ]; then
    echo "OK: ions.mdp copied to $gmxdir."
else
    echo "ERROR: ions.mdp could not be copied to $gmxdir."
    exit 1
fi
 
# Add ions with GROMACS to neutralize and reach target concentration
echo "Adding ions to neutralize the system and reach the target concentration."
gmx grompp -f "$gmxdir/ions.mdp" -c "$gmxdir/updated_gromacs.gro" -p "$gmxdir/topol.top" -o "$gmxdir/ions.tpr" -quiet -maxwarn 1
# There will be a warning because the system is not neutralized.
gmx genion -s "$gmxdir/ions.tpr" -o "$gmxdir/ions.gro" -p "$gmxdir/topol.top" -pname Na+ -nname Cl- -neutral -conc 0.15 -quiet > $gmxdir/genion.log 2>&1 << EOL
Water
quit
EOL

# Print last 3 lines of genion log
tail -n 3 $gmxdir/genion.log
echo "Ions in topology file after genion:"
grep -E "^(Na|Cl)\s+[0-9]+" "$gmxdir/topol.top"
# check if topol.top contains more than 1 Na and Cl ions
num_Na=$(grep -E "^(Na)\s+[0-9]+" "$gmxdir/topol.top" | awk '{print $2}')
num_Cl=$(grep -E "^(Cl)\s+[0-9]+" "$gmxdir/topol.top" | awk '{print $2}')
if ! [[ "$num_Na" =~ ^[0-9]+$ ]] || (( num_Na <= 1 )); then
    echo "WARNING: The number of Na ions added to the system is smaller or equal to 1. The system is likely not neutralized."
else
    echo "OK: Ions added successfully. Na: $num_Na"
fi

if ! [[ "$num_Cl" =~ ^[0-9]+$ ]] || (( num_Cl <= 1 )); then
    echo "WARNING: The number of Cl ions added to the system is smaller or equal to 1. The system is likely not neutralized."
else
    echo "OK: Ions added successfully. Cl: $num_Cl"
fi