#!/bin/bash

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
    error "gmxdir argument is required."
    print_usage
    exit 1
fi

gmxdir=$1

# check valid directory and files
check() {
    local msg="${@: -1}"
    local cmd=( "${@:1:$#-1}" )
    if ! "${cmd[@]}"; then
        error "$msg"
        exit 1
    fi
}

check [ -d "$gmxdir" ] "gmxdir '$gmxdir' does not exist or is not a directory."
check [ -f "$gmxdir/gromacs.gro" ] "gromacs.gro file not found in '$gmxdir'."
check [ -f "$gmxdir/topol.top" ] "topol.top file not found in '$gmxdir'."
check [ -d "$gmxdir/toppar" ] "toppar directory not found in '$gmxdir'."
check [ "$(find "$gmxdir/toppar" -maxdepth 1 -name '*.itp' | wc -l)" -gt 0 ] "No .itp files found in '$gmxdir/toppar'."

# Check if checkpoint file from previous run exists
if ls "$gmxdir/check_04_removeamberions_"*.done 1> /dev/null 2>&1; then
    warn "Checkpoint file found. Skipping Amber ions removal."
    # check if updated_gromacs.gro exists
    if [ -f "$gmxdir/updated_gromacs.gro" ]; then   
        ok "Updated topology file found at $gmxdir/updated_gromacs.gro."
        warn "Skipping Amber ions removal."
        exit 0
    else
        warn "Checkpoint file found but updated_gromacs.gro not found at $gmxdir. The previous run might have been incomplete. Proceeding with Amber ions removal."
    fi
fi

# count lines with Na and Cl in gromacs.gro (before removal)
num_Na=$(grep -c 'Na' "${gmxdir}/gromacs.gro")
num_Cl=$(grep -c 'Cl' "${gmxdir}/gromacs.gro")
echo "Number of Na ions in gromacs.gro: $num_Na" 
echo "Number of Cl ions in gromacs.gro: $num_Cl"

# Warnings if num_Na or num_Cl is more than 1
if [ "$num_Na" -gt 1 ] ; then
    warn "Number of Na ions is expected to be one, but was $num_Na. Later addition of ions with GROMACS may not work as expected."
fi
if [ "$num_Cl" -gt 1 ] ; then
    warn "Number of Cl ions is expected to be one, but was $num_Cl. Later addition of ions with GROMACS may not work as expected."
fi

# Find first lines containing "Na" and "Cl"
line1=$(grep -n 'Na' "${gmxdir}/gromacs.gro" | head -n 1 | cut -d: -f1)
line2=$(grep -n 'Cl' "${gmxdir}/gromacs.gro" | head -n 1 | cut -d: -f1)

last_line=$(wc -l < "${gmxdir}/gromacs.gro")

# Check if both Na and Cl lines were found
if [ -z "$line1" ] || [ -z "$line2" ]; then
    error "Could not find both Na and Cl lines" >&2
    exit 1
else
    ok "Found Na at line $line1 and Cl at line $line2 in gromacs.gro"
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

# check updated_gromacs.gro was created
if [ -f "$gmxdir/updated_gromacs.gro" ]; then   
    ok "Updated file saved as updated_gromacs.gro."
else
    error "updated_gromacs.gro could not be created."
    exit 1
fi

# check there are no Na or Cl updated_gromacs.gro
num_Na_updated=$(grep -c 'Na' "${gmxdir}/updated_gromacs.gro")
num_Cl_updated=$(grep -c 'Cl' "${gmxdir}/updated_gromacs.gro")
if [ "$num_Na_updated" -eq 0 ] && [ "$num_Cl_updated" -eq 0 ] ; then
    ok "Successfully removed Amber ions from gromacs.gro. No Na or Cl lines found in updated_gromacs.gro."
else
    warn "Found $num_Na_updated Na and $num_Cl_updated Cl lines in updated_gromacs.gro, but expected 0."
fi


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
    warn "Na remain in topology file $gmxdir/topol.top after removal. This could result in clashing nomenclature later."
fi
if [ "$remaining_Cl" -gt 0 ]; then
    warn "Cl remain in topology file $gmxdir/topol.top after removal. This could result in clashing nomenclature later."
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
        warn "Na+ remain in topology file $itpfile after removal. This could result in clashing nomenclature later."
    fi
    if [ "$remaining_Cl" -gt 0 ]; then
        warn "Cl- remain in topology file $itpfile after removal. This could result in clashing nomenclature later."
    fi
done

# Create checkpoint file
touch "$gmxdir/check_04_removeamberions_$(date +%Y%m%d_%H%M%S).done"

# TODO these were in the initial script, are they necessary?
# sed -i "s/Cl-/Cl/g" "$gmxdir/topol.top"
# sed -i "s/Na+/Na/g" "$gmxdir/topol.top"
sed -i "s/WAT/Water/g" "$gmxdir/topol.top"
remaining_WAT=$(grep -c '/WAT/d' "$gmxdir/topol.top")
warn "WAT remain in topology file $gmxdir/topol.top after removal. This could result in clashing nomenclature later."
