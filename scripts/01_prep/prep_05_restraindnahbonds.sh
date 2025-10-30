#!/bin/bash

print_usage() {
  cat <<'USAGE'
Usage: prep_05_restraindnahbonds.sh [OPTIONS]

Calculate distance restraints for DNA H-bonds and append them to a topology file.

Options:
  --gmxdir DIR        Directory containing GROMACS files          [default: .]
  --topol FILE        Topology file to append restraints to       [default: ./topol.top]
  --dna_a FILE        DNA chain A topology file (.itp)            [default: ./toppar/molecule_3.itp]
  --dna_b FILE        DNA chain B topology file (.itp)            [default: ./toppar/molecule_4.itp]
  --type INT          Restraint type                              [default: 1]
  --type2 INT         Secondary restraint type                    [default: 2]
  --low FLOAT         Lower bound distance (Å)                    [default: 0.0]
  --up1 FLOAT         Upper bound 1 (Å)                           [default: 0.3]
  --up2 FLOAT         Upper bound 2 (Å)                           [default: 0.4]
  --fac FLOAT         Force constant scaling factor               [default: 1.0]

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
	echo "ERROR: gmxdir '$gmxdir' is not a valid directory."
	print_usage
	exit 1
fi

# check topol file exists
if [ ! -f "$topol" ]; then
    echo "ERROR: topology file '$topol' does not exist."
    exit 1
fi

# Validate .itp files
if [ ! -f "$dna_a" ]; then
    echo "ERROR: DNA chain A topology file '$dna_a' not found."
    print_usage
    exit 1
fi
if [ ! -f "$dna_b" ]; then
    echo "ERROR: DNA chain B topology file '$dna_b' not found."
    print_usage
    exit 1
fi

# Validate integer parameters
if ! is_int "$type"; then
    echo "ERROR: type must be an integer. Got: '$type'"
    print_usage
    exit 1
fi
if ! is_int "$type2"; then
    echo "ERROR: type2 must be an integer. Got: '$type2'"
    print_usage
    exit 1
fi

# Validate numeric bounds
for vname in low up1 up2 fac; do
    val=${!vname}
    if ! is_num "$val"; then
        echo "ERROR: $vname must be a number. Got: '$val'"
        print_usage
        exit 1
    fi
done

# Check order low < up1 < up2
if ! float_lt "$low" "$up1"; then
    echo "ERROR: expected low < up1. Got low=$low up1=$up1"
    print_usage
    exit 1
fi
if ! float_lt "$up1" "$up2"; then
    echo "ERROR: expected up1 < up2. Got up1=$up1 up2=$up2"
    print_usage
    exit 1
fi

# fac should be positive
if ! awk -v f="$fac" 'BEGIN{ if(f<=0) exit 1; exit 0 }'; then
    echo "ERROR: fac must be > 0. Got: $fac"
    print_usage
    exit 1
fi

echo "Calculating distance restrains for h-bonds in $dna_a and $dna_b.
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

# Check equal lengths
l=${#seq_a[@]}
if [[ ${#seq_a[@]} == ${#seq_b[@]} ]]; then
	echo "OK: DNA sequences have the same length $l."
else 
	echo "ERROR: DNA sequences must have the same length."
	exit 1
fi

# Check valid DNA seq
if [[ "$seq_a" =~ [^ACGT] && "$seq_b" =~ [^ACGT] ]]; then
    echo "ERROR: Invalid DNA sequence. Only A, C, G, T are allowed."
    exit 1
fi

echo "OK: DNA sequences are valid."

# Check complementary seqs
valid_pairs=("AT" "TA" "CG" "GC")
for i in ${!seq_a[@]}; do
	a=${seq_a[$i]}
	b=${seq_b[-($i+1)]} 
	if [[ ! ${valid_pairs[@]} =~ $a$b ]]; then
		echo "ERROR: Invalid base pair. DNA sequences are not complementary."
		exit 1
	fi
done

echo "OK: DNA sequences are complementary."

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
	if [[ "$dna_a" == $file ]]; then
		echo "Match $dna_a, $file"
		global_a=$global
	elif [[ "$dna_b" == $file ]]; then
		echo "Match $dna_b, $file"
		global_b=$global
	fi
	n=$(awk '/\[ atoms \]/ {flag=1; next} /\[ bonds \]/ {flag=0} flag' $file | awk '/[0-9]/ {print $1}' | tail -n 1)
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
		echo "ERROR: Error building pairs. Invalid base."
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

echo "OK: DNA distance restraints added to $topol."
