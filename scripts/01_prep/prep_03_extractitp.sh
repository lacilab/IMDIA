#!/bin/bash

# Check if the top file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 gromacs_directory"
    exit 1
fi

# Input gromacs topology
gmxdir=$1
top="$gmxdir/topol.top"

# Extract molecule itps
n=$( awk '/\[ moleculetype \]/' $top | wc -l)

for ((i = 1 ; i <= n ; i++)); do
  awk '
    /\[ moleculetype \]/ || /\[ system \]/ {count++; if (count == 1) {print; next}} 
    count == 1 {print} 
    count == 2 {count = 3} 
    ' "$top" > "$gmxdir/toppar/molecule_$i".itp
  awk '
    /\[ moleculetype \]/ || /\[ system \]/ {count++; if (count == 1) { print "#include \"toppar/molecule_'$i'.itp\""; next}}
    count == 1 {next}
    count == 2 {count = 3}
    1
    ' "$top" > temp_file && mv temp_file "$top"
done

# Extract forcefield.itp
awk '/^#include/ {exit} 1' "$top" > "$gmxdir/toppar/forcefield.itp"
awk '
    /\[ defaults \]/  || /^#include/ {count++; if (count == 1) { print "#include \"toppar/forcefield.itp\""; next}}
    count == 1 {next}
    count == 2 {count = 3}
    1
    ' "$top" > temp_file && mv temp_file "$top"

echo "$n molecules have been extracted into separate .itp files. from $top"
