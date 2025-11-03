#!/bin/bash

# Check if the top file is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 gromacs_directory replicate_n"
    exit 1
fi

# Print usage/help
print_usage() {
    cat <<'USAGE'
Usage: md.sh gromacs_directory replicate_n
Run molecular dynamics simulations using GROMACS. Steps: energy minimization, NVT, NPT, and MD.

Positional arguments:
  gromacs_directory   Directory containing GROMACS input files (ions.gro, topol.top, etc.)
  replicate_n         Replicate number (used to create output directory)

Options:
  --with-dna          Complex includes DNA
  --with-phosphate    Complex includes phosphorylated residues
  --ntmpi N           Number of MPI ranks (default: auto-detected)
  --ntomp N           Number of OpenMP threads per rank (default: auto-detected) 
  --ngpu N            Number of GPUs to use (default: auto-detected)
  --gpulist N         Sequence of GPU IDs to use i.e. 0011 (default: auto-detected)
  --em FILE           Energy minimization parameter file (default: minim.mdp)
  --nvt FILE          NVT equilibration parameter file (default: nvt.mdp)
  --npt FILE          NPT equilibration parameter file (default: npt.mdp)
  --md FILE           Production MD parameter file (default: md.mdp)
  -h, --help          Show this help message and exit
USAGE
}

scriptdir=$(pwd)
# Positional arguments
gmxdir=$1
n=$2

# check gmxdir exists
if [ ! -d "$gmxdir" ]; then
    error "gmxdir '$gmxdir' does not exist or is not a directory."
    exit 1
fi

# if n is not provided, find repn directories in gmxdir and set n to the next available number
if [ -z "$n" ]; then
    existing_reps=($(ls -d "$gmxdir"/rep* 2>/dev/null | grep -o 'rep[0-9]*' | sed 's/rep//' | sort -n))
    if [ ${#existing_reps[@]} -eq 0 ]; then
        n=1
    else
        last_rep=${existing_reps[-1]}
        n=$((last_rep + 1))
    fi
fi

# Default options
with_dna=false
with_phosphate=false

# Resource detection. By default, use all available resources
if command -v nvidia-smi &> /dev/null; then
    NGPUS=$(nvidia-smi -L | wc -l)
else
    NGPUS=0
fi

if [ -n "$SLURM_NTASKS" ] && [ -n "$SLURM_CPUS_PER_TASK" ]; then
    NTMPI=$SLURM_NTASKS
    NTOMP=$SLURM_CPUS_PER_TASK
else
    # fallback for local machine
    NTOMP=$(nproc)
    NTMPI=1
fi

# Parse flags before positional arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-dna)
      with_dna=true
      echo "Including --with-dna flag. System should include DNA."
      shift
      ;;
    --with-phosphate)
      with_phosphate=true
      echo "Including --with-phosphate flag. System should include phosphorylated residues."
      shift
      ;;
    --ntmpi)
      NTMPI="$2"
      shift 2
      ;;
    --ntomp)
      NTOMP="$2"
      shift 2
      ;;
    --ngpu)
      NGPUS="$2"
      shift 2
      ;;
    --gpulist)
      GPULIST="$2"
      shift 2
      ;;
    --em)
      EM_FILE="$2"
      shift 2
      ;;
    --nvt)
      NVT_FILE="$2"
      shift 2
      ;;
    --npt)
      NPT_FILE="$2"
      shift 2
      ;;
    --md)
      MD_FILE="$2"
      shift 2
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

echo "Running $NTMPI MPI ranks × $NTOMP threads per rank, with $NGPUS GPUs."

if [ "$NGPUS" -gt 0 ]; then
    flag_nb_gpu="-nb gpu"
    flag_pme_gpu="-pme gpu -npme 1"
    if [ -n "$GPULIST" ]; then
        echo "Using specified GPU list: $GPULIST"
    else
        echo "No GPU list specified. Generating GPU task list automatically."
        for ((i=0;i<NTMPI;i++)); do
          gpu=$(( i % NGPUS ))   # cycle GPUs if NTMPI > NGPUS
          GPULIST+="$gpu"
        done
    fi
    flag_gpu_tasks="-gputasks $GPULIST"
else
    GPULIST=""
    echo "No GPUs detected or specified. Running on CPU only."
fi


echo "Using $NGPUS GPU(s) with gputasks: $GPULIST"

# create output dir and copy parameter files
parent_dir="$(dirname "$gmxdir")"
outdir="$parent_dir/rep${n}/gromacs"
mkdir -p "$outdir"
cp "$scriptdir"/*.mdp "$outdir"
cd "$outdir"

# check outdir exists and contains minim.mdp, nvt.mdp, npt.mdp, md.mdp
if [ ! -d "$outdir" ] || [ ! -f "minim.mdp" ] || [ ! -f "nvt.mdp" ] || [ ! -f "npt.mdp" ] || [ ! -f "md.mdp" ]; then
    error "Output directory '$outdir' is does not exist or is missing parameter files. Required files:  minim.mdp, nvt.mdp, npt.mdp, md.mdp"
    exit 1
fi

# if -em, --nvt, --npt, --md options are provided, copy those files instead of the default ones to $outdir. Else, copy the default files
# helper to set parameter files. args: <custom_path_var> <filename>
param_mdp() {
    custom="$1"
    filename="$2"
    # prefer custom if provided and exists, otherwise use default from scriptdir
    if [ -n "$custom" ]; then
        if [ -f "$custom" ]; then
            echo "Using custom $filename: $custom"
            cp -- "$custom" "$outdir/$filename"
        else
            error "Custom '$custom' not found. Using default $filename."
            cp -- "$scriptdir/$filename" "$outdir/$filename"
        fi
    else
        cp -- "$scriptdir/$filename" "$outdir/$filename"
    fi
}

# parameter files
param_mdp "$EM_FILE"  "minim.mdp"
param_mdp "$NVT_FILE" "nvt.mdp"
param_mdp "$NPT_FILE" "npt.mdp"
param_mdp "$MD_FILE"  "md.mdp"


# Energy minimization
gmx grompp -f minim.mdp -c "$gmxdir/ions.gro" -r "$gmxdir/ions.gro" -p "$gmxdir/topol.top" -o em.tpr -quiet -maxwarn 1
gmx mdrun -v -deffnm em -quiet $flag_nb_gpu

# Create groups
if [[ $with_dna == true ]] && [[ $with_phosphate == true ]]; then
gmx make_ndx -f em.gro -o index.ndx -quiet << EOL
r DA | r DC | r DT | r DG | r DC3 | r DC5 | r DG3 | r DG5
name 19 DNA
1 | 19 | 13
name 20 Protein_DNA_PTR
14 | 15 | 16
name 21 Water_and_ions 
1 | 13
name 22 Protein_PTR
quit
EOL
else
gmx make_ndx -f em.gro -o index.ndx -quiet << EOL
12 | 16 | 17
name 18 Water_and_ions 
quit
EOL
fi

# NVT
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p "$gmxdir/topol.top" -o nvt.tpr -n index.ndx -quiet -maxwarn 1
gmx_mpi mdrun -v -deffnm nvt -quiet $flag_nb_gpu

# NPT
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p "$gmxdir/topol.top" -o npt.tpr -n index.ndx -quiet -maxwarn 1 
gmx_mpi mdrun -v -deffnm npt -quiet $flag_nb_gpu

# MD
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p "$gmxdir/topol.top" -o md.tpr -n index.ndx -quiet -maxwarn 1
gmx_mpi mdrun -ntmpi $NTMPI -ntomp $NTOMP $flag_nb_gpu $flag_pme_gpu $flag_gpu_tasks -v -deffnm md -quiet > mdtrack.log 2>&1

# for restarting from a checkpoint
# gmx mdrun -ntmpi 16 -ntomp 8 -cpi md.cpt -append -nb gpu -pme gpu -npme 1 -gputasks 0000000000111111 -v -deffnm md -quiet > output2.log 2>&1 &
