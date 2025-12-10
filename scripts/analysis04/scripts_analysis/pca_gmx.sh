#!/bin/bash

# TODO remove if not needed
# source activate base
# conda activate thesis_md
# source /opt/gromacs/gromacs_2024/bin/GMXRC


# Print usage/help
print_usage() {
    cat <<'USAGE'
Usage: pca_gmx.sh --outdir OUTDIR --tpr TPRFILE --xtc XTCFILE
Run principal component analysis (PCA) on MD results using GROMACS.

Options:
  --outdir OUTDIR   Output directory for GROMACS files   (required)
  --tpr TPRFILE     TPR file for the system               (required)
  --xtc XTCFILE     XTC file for the system               (required)
  -h, --help        Show this help message and exit
USAGE
}

# if there are less than 3 arguments, print usage and exit
if [ "$#" -lt 3 ]; then
    echo "Insufficient arguments provided. Please provide --outdir, --tpr, and --xtc."
    print_usage
    exit 1
fi

# parse flags
outdir=""
gmxdir=""
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

# cd $outdir

#STEP 1: build and diagonalize the covariance matrix
gmx covar -f $xtc -s $tpr -o $outdir/eigenval.xvg -v $outdir/eigenvec.trr -av $outdir/average.pdb -l $outdir/covar.log << EOL
C-alpha
C-alpha
EOL

#STEP 2: project an MD trajectory along a 1st eigenvector
gmx anaeig -f $xtc -v $outdir/eigenvec.trr -eig $outdir/eigenval.xvg -s $tpr -proj $outdir/proj_1_2.xvg -b 0 -e 500001 -first 1 -last 1 << EOL
C-alpha
C-alpha
EOL

#STEP 3: extract time scale from proj_1_2.xvg
awk '{print $1}' $outdir/proj_1_2.xvg > $outdir/timescale.csv
vim $outdir/timescale.csv << EOL
:set number
:1,24d
:wq
EOL

#STEP 4: project an MD trajectory in two dimensions of selected eigenvectors (essential subspace)
gmx anaeig -f $xtc -v $outdir/eigenvec.trr -eig $outdir/eigenval.xvg -s $tpr -2d $outdir/2dproj_1_2.xvg  -b 0 -e 500001 -first 1 -last 2 <<EOL
C-alpha
C-alpha
EOL
vim $outdir/2dproj_1_2.xvg << EOL
:set number
:1,17d
:wq
EOL

#STEP 5: merge timescale into 2dproj_1_2.xvg
# paste 2dproj_1_2.xvg timescale.csv  > 2dproj_timescale.csv
awk 'NR==FNR {t[NR]=$1; next} {printf "%.5f,%.5f,%.5f\n", $1*10, $2*10, t[FNR]/1000}' $outdir/timescale.csv $outdir/2dproj_1_2.xvg > $outdir/2dproj_timescale.csv


#STEP 6: plot PCA 2D projection
gnuplot << EOL
set terminal png size 600, 400
set output '$outdir/2D_projection_1_2.png'
set title 'PCA 2D-projection' 
set ylabel 'PC2 (Å)'
set xlabel 'PC1 (Å)' 
set grid
unset key
set datafile separator ","
set palette defined ( 0 "#000090",\
        1 "#000fff",\
        2 "#0090ff",\
        3 "#0fffee",\
        4 "#90ff70",\
        5 "#ffee00",\
        6 "#ff7000",\
        7 "#ee0000",\
        8 "#7f0000")
set cblabel 'Time (ns)'
plot '$outdir/2dproj_timescale.csv' u 1:2:3 with points palette
quit
EOL
