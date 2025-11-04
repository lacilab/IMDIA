#!/bin/bash

source "../utils/colors.sh"

# Print usage/help
print_usage() {
    cat <<'USAGE'
Usage: prep_00_all.sh datadir filename
Prepare a system for MD.

Positional arguments:
  datadir    Directory that contains the input PDB
  filename   Name of the PDB file (with or without .pdb)

Options:
  --with-dna          Complex includes DNA
  --with-phosphate    Complex includes phosphorylated residues
  --ffs <list>        Specify force fields (space-separated list) to use in LEaP. Avoid using clashing force fields. Default: protein.ff19SB, water.opc.
  -h, --help          Show this help message and exit
USAGE
}

# Check for help argument
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    print_usage
    exit 0
fi

# Require exactly 2 positional arguments
if [ "$#" -ne 2 ]; then
    error "datadir and filename are required."
    print_usage
    exit 1
fi

########## POSITIONAL ARGUMENTS ##########

datadir="$1"
file="$2"

# Normalize filename: accept with or without .pdb, store basename without extension
file="$(basename "$file")"
if [[ "$file" == *.pdb ]]; then
    file="${file%.pdb}"
fi

# check datadir exists
if [ ! -d "$datadir" ]; then
    error "datadir '$datadir' does not exist or is not a directory."
    exit 1
fi

#### LOG

LOGFILE="${datadir}/prep_${file}_$(date +%Y%m%d_%H%M%S).log"

# Redirect all stdout and stderr to tee
# This will print to screen AND save to the logfile, removing color tags
exec > >(tee >(sed -r "s/\x1B\[[0-9;]*[mK]//g" >> "$LOGFILE")) 2>&1
echo "Logging to: $LOGFILE"

# Check the input PDB exists
if [ ! -f "$datadir/${file}.pdb" ]; then
    error "file '$datadir/${file}.pdb' not found."
    exit 1
fi

########### Check key software dependencies
# Check GROMACS is ready and print version
if ! command -v gmx >/dev/null 2>&1; then
    echo "GROMACS (gmx) not found. Please install GROMACS and ensure it is in your PATH."
    exit 1
fi
gmx --version

# Check AmberTools is ready and print version
if ! command -v ambpdb >/dev/null 2>&1; then
    echo "AmberTools (ambpdb) not found. Please install AmberTools and ensure it is in your PATH."
    exit 1
fi
ambpdb -V

# Check ParmEd is ready and print version
if ! command -v parmed >/dev/null 2>&1; then
    echo "ParmEd (parmed) not found. Please install ParmEd and ensure it is in your PATH."
    exit 1
fi
parmed --version

# Confirmation
ok "Preparing: datadir='$datadir' file='${file}.pdb'"

########## OPTIONAL ARGUMENTS ##########
# Default params
with_dna=false
with_phosphate=false
ffs=()

# Parse flags before positional arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-dna)
      with_dna=true
      echo "Including --with-dna flag for preparation. System should include DNA."
      shift
      ;;
    --with-phosphate)
      with_phosphate=true
      echo "Including --with-phosphate flag for preparation. System should include phosphorylated residues."
      shift
      ;;
    --ffs)
      shift
      # Collect all following arguments that are not options
      while [[ $# -gt 0 && ! $1 =~ ^- ]]; do
        ffs+=("$1")
        shift
      done
      echo "Custom force fields specified: ${ffs[*]}."
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

# Remove lines starting with CONECT, which are added by MOE and cause issues in AmberTools 
# check if $datadir/${file}_noconect.pdb already exists
if [ -f "$datadir/${file}_noconect.pdb" ]; then
    warn "$datadir/${file}_noconect.pdb already exists. Skipping CONECT removal."
else
    echo "Removing CONECT lines from ${file}.pdb. These lines can cause issues in AmberTools."
    sed '/^CONECT/d' "$datadir/${file}.pdb" > "$datadir/${file}_noconect.pdb"
    remaining=$(grep -c '^CONECT' "$datadir/${file}_noconect.pdb")
    echo "CONECT lines before: $(grep -c '^CONECT' "$datadir/${file}.pdb")"
    echo "CONECT lines after:  $remaining"
    if [ "$remaining" -gt 0 ]; then
        warn "$remaining CONECT lines remain in ${file}_noconect.pdb"
    fi
fi


# Check if datadir/amber already exists.
if [ -d "$datadir/amber" ]; then
    warn "$datadir/amber already exists."
    if [ -f "$datadir/amber/${file}_noconect.pdb" ]; then
        warn "$datadir/amber/${file}_noconect.pdb already exists in amber directory."
    else
        mv "$datadir/${file}_noconect.pdb" "$datadir/amber/"
    fi
else
    echo "$datadir/amber does not exist. Creating new directory."
    mkdir -p "$datadir/amber"
    mv "$datadir/${file}_noconect.pdb" "$datadir/amber/"
fi
# Check that datadir/amber exists and contains the structure
if [ -d "$datadir/amber" ] && [ -f "$datadir/amber/${file}_noconect.pdb" ]; then
    ok "$datadir/amber directory is ready. "
    ok "${file}_noconect.pdb is in $datadir/amber/ ."
else
    error "$datadir/amber directory is missing or incomplete."
    exit 1
fi

########## PREPARATION ##########

##### AmberTools preparation
# add --with-dna if $with_dna is true
((with_dna ? dna_flag="--with-dna" : dna_flag="")) 
# add --with-phosphate if $with_phosphate is true
((with_phosphate ? phosphate_flag="--with-phosphate" : phosphate_flag=""))

python prep_01_amber.py --path-pdb "${datadir}/amber/${file}_noconect.pdb"  \
    --ffs "${ffs[@]}" \
    $dna_flag \
    $phosphate_flag

# Check that datadir/gromacs_common does not exist. If it exists, show a message and backup the old one (change the name)
if [ -d "$datadir/gromacs_common" ]; then
    warn "$datadir/gromacs_common already exists."
else
    mkdir "$datadir/gromacs_common"
    gmxdir="$datadir/gromacs_common"
    if [ -d "$gmxdir" ]; then
        ok "Gromacs directory created at $gmxdir."
    else
        error "Gromacs directory could not be created at $gmxdir."
        exit 1
    fi
fi

#### Convert to Gromacs format with parmed
python prep_02_amb2gro.py -p "${datadir}/amber/${file}_noconect_protonated_solvated.prmtop" \
    -c "${datadir}/amber/${file}_noconect_protonated_solvated.inpcrd" \
    -t "${gmxdir}/topol.top" \
    -g "${gmxdir}/gromacs.gro" \


# Extract .itps from .top file
bash prep_03_extractitp.sh "$gmxdir/topol.top" $dna_flag

# Remove the 2 ions added with AmberTools from .gro files
bash prep_04_removeamberions.sh "$gmxdir"

# Add ions with GROMACS (neutralize and reach target concentration)
bash prep_06_addgmxions.sh "$gmxdir"

# Add DNA H-bonds restraints 
if [[ $with_dna == true ]]; then
    # for each .itp file, check if it contains DC, DA, DG or DT
    dna_files=()
    for file in "$gmxdir/toppar/"*.itp; do
        if grep -qE 'DC|DA|DG|DT' "$file"; then
            echo "Found DNA base in $file"
            # save $file in list of dna_files
            dna_files+=("$file")
        fi
    done
    # if dna_files does not contain 2 files, print message and exit
    if [ ${#dna_files[@]} -ne 2 ]; then
        error "--with-dna flag provided but found ${#dna_files[@]} files containing DNA residues. Expected 2 files."
        exit 1
    fi
    echo "The complex includes DNA. Adding DNA H-bond restraints."
    bash prep_05_restraindnahbonds.sh --gmxdir $gmxdir \
        --topol $gmxdir/topol.top \
        --dna_a ${dna_files[0]} \
        --dna_b ${dna_files[1]}
else    
    echo "The complex does not include DNA. Skipping DNA H-bond restraints. If this is incorrect, rerun with the --with-dna flag."
fi

ok "Preparation completed successfully."

# Create checkpoint file
touch "$datadir/check_01_prep_all_$(date +%Y%m%d_%H%M%S).done"