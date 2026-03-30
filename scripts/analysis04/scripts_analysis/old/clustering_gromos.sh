#!/bin/bash

source "../utils/colors.sh"

# Print usage/help
print_usage() {
    cat <<'USAGE'
Usage: clustering_gromos.sh --outdir OUTDIR --tpr TPRFILE --xtc XTCFILE
Run clustering using the Gromos method with GROMACS.

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

if [ ! -d "$outdir" ]; then
    mkdir "$outdir"
    echo "Creating directory: $outdir"
else
    echo "Directory already exists: $outdir"
fi

# TODO check if files exist
# TODO check that the files have the correct extensions

gmx cluster -f $xtc -s $tpr -method gromos -cutoff 0.1 \
    -o $outdir/gromos_clusters.xpm -g $outdir/gromos_clusters.log \
    -clndx $outdir/gromos_clusters.ndx -dist $outdir/gromos_clusters_rmsd_dist.xvg \
    -clid $outdir/gromos_clust_ids.xvg -sz $outdir/gromos_clust_size.xvg -quiet << EOL
C-alpha
C-alpha
EOL

gmx xpm2ps -f $outdir/gromos_clusters.xpm -o $outdir/gromos_clusters.eps -quiet

gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r300 -dEPSCrop \
   -sOutputFile=$outdir/gromos_clusters.png \
   $outdir/gromos_clusters.eps