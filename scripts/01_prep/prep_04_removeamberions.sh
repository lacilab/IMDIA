#!/bin/bash

# Check if datadir argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 gromacs_directory"
    exit 1
fi

# Print usage/help
print_usage() {
    cat <<'USAGE'
Usage: prep_04_removeamberions.sh topology
Remove Amber ions from the GROMACS topology files.

Positional arguments:
  gmxdir    Path to the GROMACS directory (containing gromacs.gro and topol.top)

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

# check that gmxdir exists and it contains gromacs.gro and topol.top and toppar/ with itp files
# TODO Less redundant
if [ ! -d "$gmxdir" ]; then
    echo "ERROR: gmxdir '$gmxdir' does not exist or is not a directory."
    exit 1
fi
if [ ! -f "$gmxdir/gromacs.gro" ]; then
    echo "ERROR: gromacs.gro file not found in '$gmxdir'."
    exit 1
fi
if [ ! -f "$gmxdir/topol.top" ]; then
    echo "ERROR: topol.top file not found in '$gmxdir'."
    exit 1
fi
if [ ! -d "$gmxdir/toppar" ]; then
    echo "ERROR: toppar directory not found in '$gmxdir'."
    exit 1
fi
if [ $(find "$gmxdir/toppar" -maxdepth 1 -name "*.itp" | wc -l) -eq 0 ]; then
    echo "ERROR: No .itp files found in '$gmxdir/toppar'."
    exit 1
fi

# First lines containing "Na" and "Cl"
line1=$(grep -n 'Na' "${gmxdir}/gromacs.gro" | head -n 1 | cut -d: -f1)
line2=$(grep -n 'Cl' "${gmxdir}/gromacs.gro" | head -n 1 | cut -d: -f1)

last_line=$(wc -l < "${gmxdir}/gromacs.gro")

# Check if both Na and Cl lines were found
if [ -z "$line1" ] || [ -z "$line2" ]; then
    echo "ERROR: Could not find both Na and Cl lines" >&2
    exit 1
fi

gawk -v last_line="$last_line" -v line1="$line1" -v line2="$line2" '
NR == 2 {atom_count = $1 - 2; print atom_count; next}  # Update the atom count on the second line
NR < line1 {                                   # Print lines before remove_line2 without changes
    print $0
    next
}
NR == line1 {next}   
NR == line2 {next}                             # Skip lines to be removed
NR >= line2 && NR <= last_line - 1 {            # Process lines after remove_line2
    # First column starts with a space, then a number followed by a string
    match($1, /([0-9]+)([A-Za-z0-9]+)/, arr);  # Match a number and the string part after it
  
    # Subtract 2 from the number
    mol_number = (arr[1] >= 2) ? arr[1] - 2 : 99998 + arr[1];
    mol_name = arr[2]
    
    # Check if the third column is a float (contains a ".")
    if ($3 !~ /\./) {
        atom_num = $3;  
        # Subtract 2, handling wrap-around
        atom_num = (atom_num >= 2) ? atom_num - 2 : 99998 + atom_num;

        printf "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", mol_number, mol_name, $2, atom_num, $4, $5, $6;
    } else {
        match($2, /(.+)(.{5})$/, arr);

        atom_str = arr[1];      # The part before the last 5 characters (string + spaces)
        atom_num = arr[2] + 0;  # Convert last 5 characters to a number (removes leading spaces)

        # Subtract 2, handling wrap-around
        atom_num = (atom_num >= 2) ? atom_num - 2 : 99998 + atom_num;

         # Print modified line
        printf "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", mol_number, mol_name, atom_str, atom_num, $3, $4, $5;
    }

}
NR == last_line {print}                               # Print the last line without changes
' "$gmxdir/gromacs.gro" > "$gmxdir/updated_gromacs.gro"


# update topology file
echo "Updating topology file to avoid nomenclature clashes in the conversion to GROMACS format."
echo "Before: Na $(grep -c '/Na/d' "$gmxdir/topol.top") and Cl $(grep -c '/Cl/d' "$gmxdir/topol.top") in $gmxdir/topol.top"
# remove Na and Cl from topol.top
sed -i "/Na/d" "$gmxdir/topol.top"
sed -i "/Cl/d" "$gmxdir/topol.top"
remaining_Na=$(grep -c '/Na/d' "$gmxdir/topol.top")
remaining_Cl=$(grep -c '/Cl/d' "$gmxdir/topol.top")
echo "After: Na $remaining_Na and Cl $remaining_Cl in $gmxdir/topol.top"
if [ "$remaining_Na" -gt 0 ]; then
    echo "WARNING: Na remain in topology file $gmxdir/topol.top after removal. This could result in clashing nomenclature later."
fi
if [ "$remaining_Cl" -gt 0 ]; then
    echo "WARNING: Cl remain in topology file $gmxdir/topol.top after removal. This could result in clashing nomenclature later."
fi

# Change ion names to gromacs standards on all .itp files in $gmxdir/toppar
for itpfile in "$gmxdir/toppar/"*.itp; do
    echo "Updating ion names in $itpfile."
    echo "Before: Na+ $(grep -c '/Na+/d' "$itpfile"), Na $(grep -c '/Na/d' "$itpfile") in $itpfile"
    echo "Before: Cl- $(grep -c '/Cl-/d' "$itpfile"), Cl $(grep -c '/Cl/d' "$itpfile") in $itpfile" 
    echo "Before: WAT $(grep -c '/WAT/d' "$itpfile"), Water $(grep -c '/Water/d' "$itpfile") in $itpfile" 
    # substitute Na+ to Na and Cl- to Cl in .itp files
    sed -i "s/Na+/Na/g" "$itpfile"
    sed -i "s/Cl-/Cl/g" "$itpfile"
    # Change the name of the solvent from WAT to Water so it can be selected as a group.
    sed -i "s/WAT/Water/g" "$itpfile"
    remaining_Na=$(grep -c '/Na+/d' "$itpfile")
    remaining_Cl=$(grep -c '/Cl-/d' "$itpfile")
    remaining_WAT=$(grep -c '/WAT/d' "$itpfile")
    echo "After: Na+ $remaining_Na, Cl- $remaining_Cl"
    if [ "$remaining_Na" -gt 0 ]; then
        echo "WARNING: Na+ remain in topology file $itpfile after removal. This could result in clashing nomenclature later."
    fi
    if [ "$remaining_Cl" -gt 0 ]; then
        echo "WARNING: Cl- remain in topology file $itpfile after removal. This could result in clashing nomenclature later."
    fi
done