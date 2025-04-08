#!/bin/bash

source activate base
conda activate thesis_md
source /opt/gromacs/gromacs_2024/bin/GMXRC

# Check usage
if [ $# -ne 1 ]; then
    echo "Usage: $0 pca_directory"
    exit 1
fi

datadir=$1
gmxdir="../../gromacs"
xtc="md_noPBC_mol_com_cluster_fit.xtc"
tpr="md.tpr"

# Copy files and change directory
cp porcupine.pml $datadir
cp modevectors.py $datadir
cd $datadir

# STEP 7: Generate extreme1.gro with gmx anaeig
gmx anaeig -v eigenvec.trr -s $gmxdir/$tpr -f $gmxdir/$xtc -b 0 -e 200001 -first 1 -last 1 -extr extreme1.gro -nframes 2001 <<EOL
C-alpha
C-alpha
EOL

# STEP 8: visualise using modevectors.py, no GUI
pymol -cq porcupine.pml

#STEP 9: clean up
rm porcupine.pml modevectors.py
rm -f \#*\#
