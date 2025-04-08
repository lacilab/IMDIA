#!/bin/bash

# Check if datadir argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 gromacs_directory"
    exit 1
fi

gmxdir=$1

# First lines containing "Na" and "Cl"
line1=$(grep -n 'Na' "${gmxdir}/gromacs.gro" | head -n 1 | cut -d: -f1)
line2=$(grep -n 'Cl' "${gmxdir}/gromacs.gro" | head -n 1 | cut -d: -f1)

last_line=$(wc -l < "${gmxdir}/gromacs.gro")

# Check if both Na and Cl lines were found
if [ -z "$line1" ] || [ -z "$line2" ]; then
    echo "Error: Could not find both Na and Cl lines" >&2
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

