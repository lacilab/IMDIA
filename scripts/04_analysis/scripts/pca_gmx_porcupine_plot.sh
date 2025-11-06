#!/bin/bash

# TODOremove if not needed
# source activate base
# conda activate thesis_md
# source /opt/gromacs/gromacs_2024/bin/GMXRC

source "../utils/colors.sh"

# Print usage/help
print_usage() {
    cat <<'USAGE'
Usage: pca_gmx_porcupine_plot.sh --outdir OUTDIR --tpr TPRFILE --xtc XTCFILE
Visualize PCA results in a porcupine plot using GROMACS and PyMOL, generating a 3D image.

Options:
  --outdir OUTDIR   Output directory for GROMACS files   (required)
  --tpr TPRFILE     TPR file for the system               (required)
  --xtc XTCFILE     XTC file for the system               (required)
  -h, --help        Show this help message and exit
USAGE
}

# if there are less than 3 arguments, print usage and exit
if [ "$#" -lt 3 ]; then
    error "Insufficient arguments provided. Please provide --outdir, --tpr, and --xtc."
    print_usage
    exit 1
fi

# parse flags
outdir=""
tpr=""
xtc=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) print_usage; exit 0 ;;
        --outdir) outdir="$2"; shift 2 ;;
        --tpr) tpr="$2"; shift 2 ;;
        --xtc) xtc="$2"; shift 2 ;; 
    esac
done

# TODO check if directories and files exist
# TODO check that the files have the correct extensions

# Copy files and change directory
cp porcupine.pml $outdir
cp modevectors.py $outdir
cd $outdir

# STEP 7: Generate extreme1.gro with gmx anaeig
gmx anaeig -v eigenvec.trr -s $tpr -f $xtc -b 0 -e 200001 -first 1 -last 1 -extr extreme1.gro -nframes 2001 <<EOL
C-alpha
C-alpha
EOL

# STEP 8: visualise using modevectors.py, no GUI
pymol -cq porcupine.pml

#STEP 9: clean up
rm porcupine.pml modevectors.py
rm -f \#*\#
