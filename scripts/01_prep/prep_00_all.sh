#!/bin/bash

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
    echo "ERROR: datadir and filename are required."
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

# Check that datadir and file exist 
if [ ! -d "$datadir" ]; then
    echo "ERROR: datadir '$datadir' does not exist or is not a directory."
    exit 1
fi

if [ ! -f "$datadir/${file}.pdb" ]; then
    echo "ERROR: file '$datadir/${file}.pdb' not found."
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
echo "OK: Preparing: datadir='$datadir' file='${file}.pdb'"

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
echo "Removing CONECT lines from ${file}.pdb. These lines can cause issues in AmberTools."
sed '/^CONECT/d' "$datadir/${file}.pdb" > "$datadir/${file}_noconect.pdb"
remaining=$(grep -c '^CONECT' "$datadir/${file}_noconect.pdb")
echo "CONECT lines before: $(grep -c '^CONECT' "$datadir/${file}.pdb")"
echo "CONECT lines after:  $remaining"
if [ "$remaining" -gt 0 ]; then
    echo "WARNING: $remaining CONECT lines remain in ${file}_noconect.pdb"
fi


# Check that datadir/amber does not exist. If it exists, show a message and backup the old one (change the name)
if [ -d "$datadir/amber" ]; then
    echo "WARNING: $datadir/amber already exists. Backing up old directory."
    mv "$datadir/amber" "$datadir/amber_backup_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$datadir/amber"
mv "$datadir/${file}_noconect.pdb" "$datadir/amber/"
# Check that datadir/amber exists and contains the structure
if [ -d "$datadir/amber" ] && [ -f "$datadir/amber/${file}_noconect.pdb" ]; then
    echo "OK: $datadir/amber directory is ready. "
    echo "OK: ${file}_noconect.pdb moved to $datadir/amber/ ."
else
    echo "ERROR: $datadir/amber directory is missing or incomplete."
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
    echo "WARNING: $datadir/gromacs_common already exists. Backing up old directory."
    mv "$datadir/gromacs_common" "$datadir/gromacs_common_backup_$(date +%Y%m%d_%H%M%S)"
fi
# Create gromacs directory
mkdir "$datadir/gromacs_common"
gmxdir="$datadir/gromacs_common"
# check gmxdir exists
if [ -d "$gmxdir" ]; then
    echo "OK: Gromacs directory created at $gmxdir."
else
    echo "ERROR: Gromacs directory could not be created at $gmxdir."
    exit 1
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
        echo "ERROR: --with-dna flag provided but found ${#dna_files[@]} files containing DNA residues. Expected 2 files."
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

echo "OK: Preparation completed successfully."