#!/bin/bash


# Print usage/help
print_usage() {
    cat <<'USAGE'
Usage: prep_03_extractitp.sh topology
Generate .itp files from topology file.

Positional arguments:
  topology    Path to the GROMACS topology file (e.g. 'topol.top')

Options:
  --with-dna          Complex includes DNA
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
    error "topology argument is required."
    print_usage
    exit 1
fi

# Input gromacs topology
top="$1"
$gmxdir=$(dirname "$top")

########## OPTIONAL ARGUMENTS ##########
# Default params
with_dna=false

# Parse flags before positional arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-dna)
      with_dna=true
      ((with_dna ? expected_itps=8 : expected_itps=6))
      echo "Including --with-dna flag for topology extraction. System should include DNA."
      shift
      ;;
    -h|--help)
      print_usage
      exit 0
      ;;
    -*)
      echo "Unknown option: $1"
      print_usage
      exit 1
      ;;
    *)
      break
      ;;
  esac
done


# Check if checkpoint file from previous run exists
if ls "$gmxdir/toppar/check_03_extractitp_"*.done 1> /dev/null 2>&1; then
    warn "Checkpoint file found. Checking .itp files."
    itps=$(find "$gmxdir/toppar" -maxdepth 1 -name "*.itp" | wc -l)
    if [ "$itps" -ne "$expected_itps" ]; then
        warn "Expected $expected_itps .itp files, but found $itps. The previous run might have been incomplete. Proceeding with .itp extraction."
    else
        ok "The number of .itp files found matches the expected count of ($expected_itps)."
        warn "Skipping .itp extraction."
        exit 0
    fi
fi

# Check that gmxdir/toppar does not exist. If it exists, show a message and backup the old one (change the name)
if [ -d "$gmxdir/toppar" ]; then
    warn "$gmxdir/toppar already exists"
else
    mkdir "$gmxdir/toppar"
    if [ -d "$gmxdir/toppar" ]; then
        ok "$gmxdir/toppar directory created at $gmxdir."
    else
        error "$gmxdir/toppar directory could not be created at $gmxdir."
        exit 1
    fi
fi

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

echo "The following .itp files have been created in $gmxdir/toppar/ and included in $top:"
$(echo $gmxdir/toppar/*.itp)

# Check that the expected number of .itp files have been created
created_itps=$(find "$gmxdir/toppar" -maxdepth 1 -name "*.itp" | wc -l)
if [ "$created_itps" -ne "$expected_itps" ]; then
    warn "Expected $expected_itps .itp files, but found $created_itps."
else
    ok "All expected .itp files have been created successfully."
    # Create checkpoint file
    touch "$gmxdir/toppar/check_03_extractitp_$(date +%Y%m%d_%H%M%S).done"
fi