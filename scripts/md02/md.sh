#!/bin/bash

source "utils/colors.sh"

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


# if there are less than two arguments, print usage and exit
if [ "$#" -lt 2 ]; then
    error "Insufficient arguments provided. gromacs_directory and replicate_n are required."
    print_usage
    exit 1
fi

##### HELPER FUNCTIONS #####

# helper to set parameter files. args: <custom_path_var> <filename>
# if -em, --nvt, --npt, --md options are provided, copy those files instead of the default ones to $outdir. Else, copy the default files
param_mdp() {
    custom="$1"
    filename="$2"
    # prefer custom if provided and exists, otherwise use default
    if [ -n "$custom" ]; then
        if [ -f "$custom" ]; then
            echo "Using custom $filename: $custom"
            cp -- "$custom" "$outdir/$filename"
        else
            error "Custom '$custom' not found. Using default $filename."
            cp -- "md02/$filename" "$outdir/$filename"
        fi
    else
        cp -- "md02/$filename" "$outdir/$filename"
    fi
}

# helper to check the step outputs. args: <step_short> <step_name>
check_step(){
  step_short="$1"
  step_name="$2"

  if grep -q "Finished mdrun" "$outdir/$step_short.log" 2>/dev/null && [ -f "$outdir/$step_short.gro" ]; then
      ok "$step_name completed successfully."
  else
      error "$step_name failed. Check $outdir/$step_short.log for details."
      exit 1
  fi
}


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

# create output dir
parent_dir="$(dirname "$gmxdir")"
outdir="$parent_dir/rep${n}/gromacs"
# check if outdir exists, if not create it
if [ ! -d "$outdir" ]; then
    echo "Creating output directory: $outdir"
    mkdir -p "$outdir"
fi  

# Check for checkpoint files
if [ -f "${outdir}/check_02_md_*.done" ]; then
    warn "Checkpoint file found. Skipping MD steps: minimization, NVT, NPT and production MD."
    echo "To rerun these steps, please remove the checkpoint file(s): ${outdir}/check_02_md_*.done"
    exit 0
fi

#### LOG
LOGFILE="${outdir}/mdsteps_$rep{n}_$(date +%Y%m%d_%H%M%S).log"
# Redirect all stdout and stderr to tee
# This will print to screen AND save to the logfile, removing color tags
exec > >(tee >(sed -r "s/\x1B\[[0-9;]*[mK]//g" >> "$LOGFILE")) 2>&1
echo "Logging to: $LOGFILE"


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

# Parse flags 
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

##### RESOURCE SETUP #####
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

##### PARAMETER FILES SETUP #####
# copy parameter files to outdir

# parameter files
param_mdp "$EM_FILE"  "minim.mdp"
param_mdp "$NVT_FILE" "nvt.mdp"
param_mdp "$NPT_FILE" "npt.mdp"
param_mdp "$MD_FILE"  "md.mdp"

# check outdir exists and contains minim.mdp, nvt.mdp, npt.mdp, md.mdp
if [ ! -d "$outdir" ] || [ ! -f "$outdir/minim.mdp" ] || [ ! -f "$outdir/nvt.mdp" ] || [ ! -f "$outdir/npt.mdp" ] || [ ! -f "$outdir/md.mdp" ]; then
    error "Output directory '$outdir' does not exist or is missing parameter files. Required files:  minim.mdp, nvt.mdp, npt.mdp, md.mdp"
    exit 1
fi

echo "Working directory: $(pwd), gmxdir: $gmxdir, outdir: $outdir"

###### Energy minimization
# if em.log exists and contains "Finished mdrun", and em.gro exists, skip minimization
if grep -q "Finished mdrun" "$outdir/em.log" 2>/dev/null && [ -f "$outdir/em.gro" ]; then
  warn "em.log indicates minimization already completed and em.gro exists. Skipping energy minimization."
else
  echo "Starting energy minimization."
  gmx grompp -f "$outdir/minim.mdp" -c "$gmxdir/ions.gro" -r "$gmxdir/ions.gro" -p "$gmxdir/topol.top" -o "$outdir/em.tpr" -quiet -maxwarn 1
  gmx mdrun -v -deffnm "$outdir/em" -quiet $flag_nb_gpu
fi

check_step "em" "Energy minimization"

# Create groups
if [[ $with_dna == true ]] && [[ $with_phosphate == true ]]; then
  echo "DNA and PHOSPHATE present: creating combined groups."
  gmx make_ndx -f "$outdir/em.gro" -o "$outdir/index.ndx" -quiet << EOL
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
elif [[ $with_dna == true ]]; then
  echo "DNA present: creating combined groups."
gmx make_ndx -f "$outdir/em.gro" -o "$outdir/index.ndx" -quiet << EOL
1 | 12
name 16 Protein_DNA
13 | 14 | 15
name 17 Water_and_ions 
quit
EOL
else
gmx make_ndx -f "$outdir/em.gro" -o "$outdir/index.ndx" -quiet << EOL
12 | 16 | 17
name 18 Water_and_ions 
quit
EOL
fi


###### NVT Equilibration
if grep -q "Finished mdrun" "$outdir/nvt.log" 2>/dev/null && [ -f "$outdir/nvt.gro" ]; then
  warn "nvt.log indicates NVT equilibration already completed and nvt.gro exists. Skipping NVT equilibration."
else
  echo "Starting NVT equilibration."
  gmx grompp -f "$outdir/nvt.mdp" -c "$outdir/em.gro" -r "$outdir/em.gro" -p "$gmxdir/topol.top" -o "$outdir/nvt.tpr" -n "$outdir/index.ndx" -quiet -maxwarn 1
  gmx_mpi mdrun -v -deffnm "$outdir/nvt" -quiet $flag_nb_gpu
fi

check_step "nvt" "NVT equilibration"


###### NPT Equilibration
if grep -q "Finished mdrun" "$outdir/npt.log" 2>/dev/null && [ -f "$outdir/npt.gro" ]; then
  warn "npt.log indicates NPT equilibration already completed and npt.gro exists. Skipping NPT equilibration."
else
  echo "Starting NPT equilibration."
  gmx grompp -f "$outdir/npt.mdp" -c "$outdir/nvt.gro" -r "$outdir/nvt.gro" -t "$outdir/nvt.cpt" -p "$gmxdir/topol.top" -o "$outdir/npt.tpr" -n "$outdir/index.ndx" -quiet -maxwarn 1 
  gmx_mpi mdrun -v -deffnm "$outdir/npt" -quiet $flag_nb_gpu
fi

check_step "npt" "NPT equilibration"

###### Production MD
if grep -q "Finished mdrun" "$outdir/md.log" 2>/dev/null && [ -f "$outdir/md.gro" ]; then
  warn "md.log indicates production MD already completed and md.gro exists. Skipping production MD."
elif [ -f "$outdir/md.cpt" ]; then
  echo "Restarting production MD from checkpoint md.cpt."
  gmx_mpi mdrun -ntmpi $NTMPI -ntomp $NTOMP $flag_nb_gpu $flag_pme_gpu $flag_gpu_tasks \
    -cpi "$outdir/md.cpt" -append -v -deffnm "$outdir/md" -quiet > "$outdir/mdtrack_$(date +%Y%m%d_%H%M%S).log" 2>&1
else
  echo "Starting production MD."
  gmx grompp -f "$outdir/md.mdp" -c "$outdir/npt.gro" -r "$outdir/npt.gro" -t "$outdir/npt.cpt" -p "$gmxdir/topol.top" -o "$outdir/md.tpr" -n "$outdir/index.ndx" -quiet -maxwarn 1
  gmx_mpi mdrun -ntmpi $NTMPI -ntomp $NTOMP $flag_nb_gpu $flag_pme_gpu $flag_gpu_tasks \
    -v -deffnm "$outdir/md" -quiet > "$outdir/mdtrack_$(date +%Y%m%d_%H%M%S).log" 2>&1
fi

check_step "md" "Production MD"

# Create checkpoint file
touch "${outdir}/check_02_md_$(date +%Y%m%d_%H%M%S).done"