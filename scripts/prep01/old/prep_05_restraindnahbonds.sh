#!/bin/bash

print_usage() {
  cat <<'USAGE'
Usage: prep_05_restraindnahbonds.sh [OPTIONS]

Calculate distance restraints for DNA H-bonds and append them to a topology file.

Options:
  --gmxdir DIR          Directory containing GROMACS files          [default: .]
  --topol FILE          Topology file to append restraints to       [default: ./topol.top]
  --dna_a FILE          DNA chain A topology file (.itp)            [default: ./toppar/molecule_3.itp]
  --dna_b FILE          DNA chain B topology file (.itp)            [default: ./toppar/molecule_4.itp]
  --range_a INT,INT     Range of residues to consider for DNA A [default: full length]
  --range_b INT,INT     Range of residues to consider for DNA B [default: full length]
  --type INT            Restraint type                              [default: 1]
  --type2 INT           Secondary restraint type                    [default: 2]
  --low FLOAT           Lower bound distance (Å)                    [default: 0.0]
  --up1 FLOAT           Upper bound 1 (Å)                           [default: 0.3]
  --up2 FLOAT           Upper bound 2 (Å)                           [default: 0.4]
  --fac FLOAT           Force constant scaling factor               [default: 1.0]

General:
  -h, --help          Show this help message and exit

USAGE
}

# PARAMETERS - Default values
# Files
gmxdir="."
topol="./topol.top"
dna_a="./toppar/molecule_3.itp"
dna_b="./toppar/molecule_4.itp"
range_a_start=1
range_a_end=-1
range_b_start=1
range_b_end=-1
# Restraints
type=1
type2=2
low=0.0
up1=0.3
up2=0.4
fac=1.0

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --gmxdir) gmxdir="$2"; shift 2 ;;
    --topol) topol="$2"; shift 2 ;;
    --dna_a) dna_a="$2"; shift 2 ;;
    --dna_b) dna_b="$2"; shift 2 ;;
    --range_a) range_a="$2"; IFS=',' read -r range_a_start range_a_end <<< "$range_a"; echo "Range A: $range_a"; shift 2 ;;
    --range_b) range_b="$2"; IFS=',' read -r range_b_start range_b_end <<< "$range_b"; echo "Range B: $range_b"; shift 2 ;;
    --type) type="$2"; shift 2 ;;
    --type2) type2="$2"; shift 2 ;;
    --low) low="$2"; shift 2 ;;
    --up1) up1="$2"; shift 2 ;;
    --up2) up2="$2"; shift 2 ;;
    --fac) fac="$2"; shift 2 ;;
    -h|--help)
      print_usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      print_usage
      exit 1
      ;;
  esac
done


# quick -h/--help
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    print_usage
    exit 0
fi

# Validate gmxdir
if [ ! -d "$gmxdir" ]; then	
	error "gmxdir '$gmxdir' is not a valid directory."
	print_usage
	exit 1
fi

# check topol file exists
if [ ! -f "$topol" ]; then
    error "topology file '$topol' does not exist."
    exit 1
fi

# Validate .itp files
if [ ! -f "$dna_a" ]; then
    error "DNA chain A topology file '$dna_a' not found."
    print_usage
    exit 1
fi
if [ ! -f "$dna_b" ]; then
    error "DNA chain B topology file '$dna_b' not found."
    print_usage
    exit 1
fi

is_int() { [[ $1 =~ ^[0-9]+$ ]]; }

is_float() {
    [[ $1 =~ ^-?[0-9]+([.][0-9]+)?$ ]]; 
}


# Validate integer parameters
if ! is_int "$type"; then
    error "type must be an integer. Got: '$type'"
    print_usage
    exit 1
fi
if ! is_int "$type2"; then
    error "type2 must be an integer. Got: '$type2'"
    print_usage
    exit 1
fi

# Validate numeric bounds
for vname in low up1 up2 fac; do
    val=${!vname}
    if ! is_float "$val"; then
        error "$vname must be a number. Got: '$val'"
        print_usage
        exit 1
    fi
done

# Check order low < up1 < up2
if (( $(echo "$low > $up1" | bc -l) )); then
    error "expected low < up1. Got low=$low up1=$up1"
    print_usage
    exit 1
fi

if (( $(echo "$up1 > $up2" | bc -l) )); then
    error "expected up1 < up2. Got up1=$up1 up2=$up2"
    print_usage
    exit 1
fi


# fac should be positive
if ! awk -v f="$fac" 'BEGIN{ if(f<=0) exit 1; exit 0 }'; then
    error "fac must be > 0. Got: $fac"
    print_usage
    exit 1
fi

# Check if checkpoint file from previous run exists
if ls "$gmxdir/check_05_restraindnahbonds_"*.done 1> /dev/null 2>&1; then
    warn "Checkpoint file found. Skipping DNA H-bond restraints addition."
    # Check if topol file contains distance_restraints section
    if grep -q "\[ distance_restraints \]" "$topol"; then
        ok "[ distance_restraints ] section found in $topol."
        warn "Skipping DNA H-bond restraints addition."
        exit 0
    else
        warn "Checkpoint file found but [ distance_restraints ] section not found in $topol. The previous run might have been incomplete. Proceeding with DNA H-bond restraints addition."
    fi
fi

# Start calculating distance restraints
echo "Calculating distance restraints for h-bonds in $dna_a and $dna_b.
Restrains will be writen to $topol.
Parameters: 
type    index    type'    low    up1    up2    fac
$type	i	 $type2	  $low	 $up1	$up2   $fac"


# Extract lines between [atoms] and [bonds]
awk '/\[ atoms \]/ {flag=1; next} /\[ bonds \]/ {flag=0} flag' $dna_a | grep -v ^";" > dna_a_atoms.txt 
readarray -t seq_a < <(awk '!seen[$3]++ {print substr($4, 2, 1)}' dna_a_atoms.txt)
# remove empty elements
seq_a=($(echo "${seq_a[@]}" | tr -s '\n' ' ' | tr -s ' ')) 
echo "Sequence A: ${seq_a[@]}"

awk '/\[ atoms \]/ {flag=1; next} /\[ bonds \]/ {flag=0} flag' $dna_b | grep -v ^";" > dna_b_atoms.txt 
readarray -t seq_b < <(awk '!seen[$3]++ {print substr($4, 2, 1)}' dna_b_atoms.txt)
# remove empty elements
seq_b=($(echo "${seq_b[@]}" | tr -s '\n' ' ' | tr -s ' ')) 
echo "Sequence B: ${seq_b[@]}"


# Apply residue range filters if provided
if [[ $range_a_start -gt 1 || $range_a_end -ne -1 ]]; then
    echo "Applying residue range filter for DNA A: $range_a_start to $range_a_end"
    if [[ $range_a_end -eq -1 ]]; then
        range_a_end=${#seq_a[@]}
    fi
    seq_a=("${seq_a[@]:$((range_a_start-1)):$((range_a_end-range_a_start+1))}")
    echo "Filtered Sequence A: ${seq_a[@]}"
fi

if [[ $range_b_start -gt 1 || $range_b_end -ne -1 ]]; then
    echo "Applying residue range filter for DNA B: $range_b_start to $range_b_end"
    if [[ $range_b_end -eq -1 ]]; then
        range_b_end=${#seq_b[@]}
    fi
    seq_b=("${seq_b[@]:$((range_b_start-1)):$((range_b_end-range_b_start+1))}")
    echo "Filtered Sequence B: ${seq_b[@]}"
fi

# Check equal lengths
l=${#seq_a[@]}
if [[ ${#seq_a[@]} == ${#seq_b[@]} ]]; then
	ok "DNA sequences have the same length $l."
else 
	error "DNA sequences must have the same length."
	exit 1
fi

# Check valid DNA seq
if [[ "$seq_a" =~ [^ACGT] && "$seq_b" =~ [^ACGT] ]]; then
    error "Invalid DNA sequence. Only A, C, G, T are allowed."
    exit 1
fi

ok "DNA sequences are valid."

# Check complementary seqs
valid_pairs=("AT" "TA" "CG" "GC")
for i in ${!seq_a[@]}; do
	a=${seq_a[$i]}
	b=${seq_b[-($i+1)]} 
	if [[ ! ${valid_pairs[@]} =~ $a$b ]]; then
		error "Invalid base pair. DNA sequences are not complementary."
		exit 1
	fi
done

ok "DNA sequences are complementary."

# Add headers of the restrain section to topol file
echo "
[ intermolecular_interactions ]
[ distance_restraints ]
; ai    aj    type    index    type'    low    up1    up2    fac" >> $topol

# Initialize dna_retraints_pairs.txt
echo "## Atom pairs" > dna_retraints_pairs.txt

# Local ids to global id
files=$(awk '/^#include/ && !/forcefield/ {print $2}' $topol | tr -d '"')
echo "Found .itp files included in $topol (excluding forcefield.itp):"
echo $files
global=0
for file in $files; do
	if [[ "$dna_a" == "$gmxdir/$file" ]]; then
		echo "Match $dna_a, $file"
		global_a=$global
	elif [[ "$dna_b" == "$gmxdir/$file" ]]; then
		echo "Match $dna_b, $file"
		global_b=$global
	fi
	n=$(awk '/\[ atoms \]/ {flag=1; next} /\[ bonds \]/ {flag=0} flag' "$gmxdir/$file" | awk '/[0-9]/ {print $1}' | tail -n 1)
	echo "$n atoms in $file"
	global=$(($global + $n))
done

echo "Number of atoms in the complex before DNA A: $global_a"
echo "Number of atoms in the complex before DNA B: $global_b"

for i in ${!seq_a[@]}; do
	a=${seq_a[$i]}
	b=${seq_b[-($i+1)]}
	case "$a" in
	    "A")
			atoms_a=("N6" "N1")
			atoms_b=("O4" "N3")
			comp="DT"
		;;
	    "T")
			atoms_a=("O4" "N3")
			atoms_b=("N6" "N1")
			comp="DA"
		;;
	    "C")
			atoms_a=("O2" "N3" "N4")
			atoms_b=("N2" "N1" "O6")
			comp="DG"
		;;
	    "G")
			atoms_a=("N2" "N1" "O6")
			atoms_b=("O2" "N3" "N4")
			comp="DC"
		;;
	    *)
		error "Error building pairs. Invalid base."
		exit 1
		;;
	esac
	for atom in ${!atoms_a[@]}; do
		# Search .itp for gromacs atom ids
		atom_a=$(awk -v i="$((i+1))" -v a="D$a" -v atom="${atoms_a[$atom]}" '$3==i && index($4, a) && $5==atom { print $1 }' dna_a_atoms.txt)
		atom_b=$(awk -v i="$((l-i))" -v comp="$comp" -v atom="${atoms_b[$atom]}" '$3==i && index($4, comp) && $5==atom { print $1 }' dna_b_atoms.txt)
		# Add atom pair to distance restrains matrix
		echo -e "$((atom_a+global_a))\t$((atom_b+global_b))\t$type\t$((i+1))\t$type2\t$low\t$up1\t$up2\t$fac" >> $topol
		# Add pairs with rsnr and atom pairs as "DNn atom DNn atom"
		echo -e "D${a^^}$((i+1)) ${atoms_a[$atom]} \t D${b^^}$((l-i)) ${atoms_b[$atom]}" >> dna_retraints_pairs.txt
	done
	
done
rm dna_a_atoms.txt dna_b_atoms.txt dna_retraints_pairs.txt

ok "DNA distance restraints added to $topol."

# Create checkpoint file
touch "$gmxdir/check_05_restraindnahbonds_$(date +%Y%m%d_%H%M%S).done"
