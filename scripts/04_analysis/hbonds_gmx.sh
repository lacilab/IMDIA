#!/bin/bash
source /opt/gromacs/gromacs_2024/bin/GMXRC

mut=$1
sample=$2

gmxdir="../../data/1bf5/$mut/$sample/gromacs"
cd "$gmxdir"

outdir="../analysis/002.hbonds"

if [ ! -d "$outdir" ]; then
    mkdir "$outdir"
    echo "Creating directory: $outdir"
else
    echo "Directory already exists: $outdir"
fi


# Create group Protein_PTR
gmx make_ndx -f md.gro -n index.ndx -o indexmd.ndx -quiet << EOL
1 | 13
quit
EOL


# Select 22-Protein_PTR and 19-DNA
gmx hbond -f md_noPBC_mol_com_cluster_fit.xtc -s md.tpr -n indexmd.ndx -b 100000 \
            -o ${outdir}/hbonds_gmx.ndx -pf yes -num ${outdir}/hbonds_gmx_num.xvg << EOL
22
19
EOL