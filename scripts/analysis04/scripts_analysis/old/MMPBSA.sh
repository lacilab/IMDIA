#!/bin/bash

# TODO remove if not needed
# source "$(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh"
# conda activate gmxMMPBSA

# Print usage/help
print_usage() {
    cat <<'USAGE'
Usage: clustering_gromos.sh --outdir OUTDIR --tpr TPRFILE --xtc XTCFILE
Run clustering using the Gromos method with GROMACS.

Options:
  --scriptdir SCRIPTDIR   Directory for scripts             (required)
  --gmxdir GMXDIR   GROMACS directory
  --outdir OUTDIR   Output directory for GROMACS files    (required)
  --tpr TPRFILE     TPR file for the system               (required)
  --xtc XTCFILE     XTC file for the system               (required)
  --top TOPFILE     TOP file for the system               (required)
  -h, --help        Show this help message and exit
USAGE
}

if [ "$#" -lt 6 ]; then
    echo "Insufficient arguments provided. Please provide --tpr, --xtc, --top, --gmxdir, --outdir, and --scriptdir."
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
        --scriptdir) scriptdir="$2"; shift 2 ;;
        --gmxdir) gmxdir="$2"; shift 2 ;;
        --outdir) outdir="$2"; shift 2 ;;
        --tpr) tpr="$2"; shift 2 ;;
        --xtc) xtc="$2"; shift 2 ;; 
        --top) top="$2"; shift 2 ;;
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

cp "$scriptdir/mmpbsa.in" "$outdir"

mpirun --use-hwthread-cpus -np 4 gmx_MMPBSA -i $outdir/mmpbsa.in -cs $tpr -ci ${gmxdir}/index.ndx -cg Chain_A Chain_B \
    -ct $xtc -cp $top -nogui -eo $outdir/FINAL_ENERGY.csv -deo $outdir/FINAL_ENERGY_DECOMP.csv

