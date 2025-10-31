#!/bin/bash

# Print usage/help
print_usage() {
    cat <<'USAGE'
Usage: prep_00_all_positional.sh datadir filename
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

# Check that gmxdir/toppar does not exist. If it exists, show a message and backup the old one (change the name)
if [ -d "$gmxdir/toppar" ]; then
    echo "WARNING: $gmxdir/toppar already exists. Backing up old directory."
    mv "$gmxdir/toppar" "$gmxdir/toppar_backup_$(date +%Y%m%d_%H%M%S)"
fi
mkdir "$gmxdir/toppar"
if [ -d "$gmxdir/toppar" ]; then
    echo "OK: $gmxdir/toppar directory created at $gmxdir."
else
    echo "ERROR: $gmxdir/toppar directory could not be created at $gmxdir."
    exit 1
fi

# Extract itps
bash prep_03_extractitp.sh "$gmxdir"

# Remove the 2 ions added with AmberTools from .gro files
bash prep_04_removeamberions.sh "$gmxdir"

# TODO move these to prep_04_removeamberions.sh?
# update topology file
echo "Updating topology file avoid nomenclature clashes in the conversion to GROMACS format."
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


# Check the ions.mdp exists
if [ ! -f ../02_md/ions.mdp ]; then
    echo "ERROR: ../02_md/ions.mdp file not found."
    exit 1
fi
cp ../02_md/ions.mdp $gmxdir
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
// check if topol.top contains Na and the number after Na is bigger than 1
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

# Add DNA H-bonds constrains
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


# Change ion names to gromacs standards
# TODO include these earlier? Are some redundant?
sed -i "s/Cl-/Cl/g" topol.top
sed -i "s/Na+/Na/g" topol.top
sed -i "s/Water/WAT/g" toppar/molecule_7.itp toppar/forcefield.itp topol.top
# TODO Number of Cl-, N+, Water changed. Check these are the correct files and the only ones. 