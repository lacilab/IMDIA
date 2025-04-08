#!/bin/bash

# Get input (seqs.itp, 5'-3' for both strands)
itp_a=${1:-toppar/molecule_3.itp} # DNAA
itp_b=${2:-toppar/molecule_4.itp} # DNAB
topol=${3:-topol.top}
type=${4:-1}
type2=${5:-2}
low=${6:-0.0}
up1=${7:-0.3}
up2=${8:-0.4}
fac=${9:-1.0}

echo "Calculating distance restrains for h-bonds in $itp_a and $itp_b.
Restrains will be writen to $topol.
Parameters: 
type    index    type'    low    up1    up2    fac
$type	i	 $type2	  $low	 $up1	$up2   $fac"

# Extract lines between [atoms] and [bonds]
awk '/\[ atoms \]/ {flag=1; next} /\[ bonds \]/ {flag=0} flag' $itp_a | grep -v ^";" > dna_a_atoms.txt 
readarray -t seq_a < <(awk '!seen[$3]++ {print substr($4, 2, 1)}' dna_a_atoms.txt)
# remove empty elements
seq_a=($(echo "${seq_a[@]}" | tr -s '\n' ' ' | tr -s ' ')) 
echo "Sequence A: ${seq_a[@]}"

awk '/\[ atoms \]/ {flag=1; next} /\[ bonds \]/ {flag=0} flag' $itp_b | grep -v ^";" > dna_b_atoms.txt 
readarray -t seq_b < <(awk '!seen[$3]++ {print substr($4, 2, 1)}' dna_b_atoms.txt)
# remove empty elements
seq_b=($(echo "${seq_b[@]}" | tr -s '\n' ' ' | tr -s ' ')) 
echo "Sequence B: ${seq_b[@]}"

# Check equal lengths
l=${#seq_a[@]}
if [[ ${#seq_a[@]} == ${#seq_b[@]} ]]; then
	echo "Seqs have the same length $l"
else 
	echo "Seqs of different lengths"
	exit 1
fi

# Check valid DNA seq
if [[ "$seq_a" =~ [^ACGT] && "$seq_b" =~ [^ACGT] ]]; then
    echo "Invalid seqs"
    exit 1
fi

echo "Seqs are valid DNA."

# Check complementary seqs
valid_pairs=("AT" "TA" "CG" "GC")
for i in ${!seq_a[@]}; do
	a=${seq_a[$i]}
	b=${seq_b[-($i+1)]} 
	if [[ ! ${valid_pairs[@]} =~ $a$b ]]; then
		echo "Invalid pair. Seqs are not complementary"
		exit 1
	fi
done

echo "Seqs are complementary."

# File headers
echo "
[ intermolecular_interactions ]
[ distance_restraints ]
; ai    aj    type    index    type'    low    up1    up2    fac" >> $topol

echo "## Atom pairs" > pairs.txt

# Local ids to global id
files=$(awk '/^#include/ && !/forcefield/ {print $2}' $topol | tr -d '"')
echo "Found .itp files included in $topol (excluding forcefield.itp):"
echo $files
global=0
for file in $files; do
	if [[ "$itp_a" == $file ]]; then
		echo "Match $itp_a, $file"
		global_a=$global
	elif [[ "$itp_b" == $file ]]; then
		echo "Match $itp_b, $file"
		global_b=$global
	fi
	n=$(awk '/\[ atoms \]/ {flag=1; next} /\[ bonds \]/ {flag=0} flag' $file | awk '/[0-9]/ {print $1}' | tail -n 1)
	echo "$n atoms in $file"
	global=$(($global + $n))
done

echo "Number of atoms before DNA A: $global_a"
echo "Number of atoms before DNA B: $global_b"

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
		echo "Error building pairs. Invalid base."
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
		echo -e "D${a^^}$((i+1)) ${atoms_a[$atom]} \t D${b^^}$((l-i)) ${atoms_b[$atom]}" >> pairs.txt
	done
	
done
rm dna_a_atoms.txt dna_b_atoms.txt pairs.txt

echo "Restrains added to $topol."
