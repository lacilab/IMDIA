#!/bin/bash

source activate base
conda activate thesis_md
source /opt/gromacs/gromacs_2024/bin/GMXRC

if [ $# -ne 1 ]; then
    echo "Usage: $0 pca_directory"
    exit 1
fi

gmxdir="../../gromacs"
xtc="md_noPBC_mol_com_cluster_fit.xtc"
tpr="md.tpr"

cd $1

#STEP 1: build and diagonalize the covariance matrix
gmx covar -f $gmxdir/$xtc -s $gmxdir/$tpr << EOL
C-alpha
C-alpha
EOL

#STEP 2: project an MD trajectory along a 1st eigenvector
gmx anaeig -f $gmxdir/$xtc -v eigenvec.trr -eig eigenval.xvg -s $gmxdir/$tpr -proj proj_1_2.xvg -b 0 -e 500001 -first 1 -last 1 << EOL
C-alpha
C-alpha
EOL

#STEP 3: extract time scale from proj_1_2.xvg
awk '{print $1}' proj_1_2.xvg > timescale.csv
vim timescale.csv << EOL
:set number
:1,24d
:wq
EOL

#STEP 4: project an MD trajectory in two dimensions of selected eigenvectors (essential subspace)
gmx anaeig -f $gmxdir/$xtc -v eigenvec.trr -eig eigenval.xvg -s $gmxdir/$tpr -2d 2dproj_1_2.xvg  -b 0 -e 500001 -first 1 -last 2 <<EOL
C-alpha
C-alpha
EOL
vim 2dproj_1_2.xvg << EOL
:set number
:1,17d
:wq
EOL

#STEP 5: merge timescale into 2dproj_1_2.xvg
# paste 2dproj_1_2.xvg timescale.csv  > 2dproj_timescale.csv
awk 'NR==FNR {t[NR]=$1; next} {printf "%.5f,%.5f,%.5f\n", $1*10, $2*10, t[FNR]/1000}' timescale.csv 2dproj_1_2.xvg > 2dproj_timescale.csv


#STEP 6: plot PCA 2D projection
gnuplot << EOL
set terminal png size 600, 400
set output '2D_projection_1_2.png'
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
plot '2dproj_timescale.csv' u 1:2:3 with points palette
quit
EOL
