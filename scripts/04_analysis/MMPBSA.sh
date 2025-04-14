#!/bin/bash

source activate base
conda activate gmxMMPBSA

if [ $# -ne 5 ]; then
    echo "Usage: $0 tpr xtc top index outdir"
    exit 1
fi

tpr=$1
xtc=$2
top=$3
index=$4
outdir=$5
oridir=$(pwd)

if [ ! -d "$outdir" ]; then
    mkdir "$outdir"
    echo "Creating directory: $outdir"
else
    echo "Directory already exists: $outdir"
fi

cp mmpbsa.in "$outdir"
cd "$outdir"

mpirun --use-hwthread-cpus -np 12 gmx_MMPBSA -i mmpbsa.in -cs "../../gromacs/$tpr" -ci "../../gromacs/$index" -cg DNA Protein_PTR \
    -ct "../../gromacs/$xtc" -cp "../../../gromacs_common/$top" -nogui -eo FINAL_ENERGY.csv -deo FINAL_ENERGY_DECOMP.csv


rm mmpbsa.in
