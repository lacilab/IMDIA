#!/bin/bash

# TODO remove if not needed
# source "$(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh"
# conda activate gmxMMPBSA

source "../utils/colors.sh"

# Print usage/help
print_usage() {
    cat <<'USAGE'
Usage: clustering_gromos.sh --outdir OUTDIR --tpr TPRFILE --xtc XTCFILE
Run clustering using the Gromos method with GROMACS.

Options:
  --tpr TPRFILE     TPR file for the system               (required)
  --xtc XTCFILE     XTC file for the system               (required)
  --top TOPFILE     TOP file for the system               (required)
  --index NDXFILE   Index file for the system             (required)
  --outdir OUTDIR   Output directory for GROMACS files    (required)
  -h, --help        Show this help message and exit
USAGE
}

if [ "$#" -lt 5 ]; then
    error "Insufficient arguments provided. Please provide --tpr, --xtc, --top, --index, and --outdir."
    print_usage
    exit 1
fi

# parse flags
outdir=""
tpr=""
xtc=""
top=""
index=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) print_usage; exit 0 ;;
        --outdir) outdir="$2"; shift 2 ;;
        --tpr) tpr="$2"; shift 2 ;;
        --xtc) xtc="$2"; shift 2 ;; 
        --top) top="$2"; shift 2 ;;
        --index) index="$2"; shift 2 ;;
        -*)
            error "Unknown option: $1"
            print_usage
            exit 1
            ;;
        *) break ;;

    esac
done


if [ ! -d "$outdir" ]; then
    mkdir "$outdir"
    echo "Creating directory: $outdir"
else
    echo "Directory already exists: $outdir"
fi

# TODO check files exist and have correct extensions

cp mmpbsa.in "$outdir"

mpirun --use-hwthread-cpus -np 12 gmx_MMPBSA -i mmpbsa.in -cs $tpr -ci $index -cg DNA Protein_PTR \
    -ct $xtc -cp $top -nogui -eo $outdir/FINAL_ENERGY.csv -deo $outdir/FINAL_ENERGY_DECOMP.csv

