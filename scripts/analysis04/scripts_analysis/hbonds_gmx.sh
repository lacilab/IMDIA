#!/bin/bash

source "../../utils/colors.sh"

# Print usage/help
print_usage() {
    cat <<'USAGE'
Usage: hbonds_gmx.sh --struct STRUCT --mut MUT --rep REP
Calculate hydrogen bonds using GROMACS.

Options:
  --struct STRUCT   Structure identifier (e.g., 1bf5)   (required)
  --mut MUT         Mutation identifier (e.g., WT)     (required)
  --rep REP         Sample/replicate identifier (e.g., rep1)     (required)
  -h, --help        Show this help message and exit
USAGE
}

# if there are less than three arguments, print usage and exit
if [ "$#" -lt 3 ]; then
    error "Insufficient arguments provided. Please provide --struct, --mut and --rep."
    print_usage
    exit 1
fi

# parse flags
struct=""
mut=""
rep=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) print_usage; exit 0 ;;
        --struct) struct="$2"; shift 2 ;;
        --mut)    mut="$2";    shift 2 ;;
        --rep) rep="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; print_usage; exit 1 ;;
    esac
done


gmxdir="data/$struct/$mut/$rep/gromacs"
cd "$gmxdir"

outdir="data/$struct/$mut/$rep/analysis/002_hbonds"

if [ ! -d "$outdir" ]; then
    mkdir "$outdir"
    ok "Creating directory: $outdir"
else
    warn "Directory already exists: $outdir"
fi

gmx make_ndx -f ${gmxdir}/md_sys.pdb -o ${gmxdir}/index_test.ndx << EOL
chain A
name 18 Chain_A
chain B
name 19 Chain_B
q
EOL

# Select groups
gmx hbond -f ${gmxdir}/md_noPBC_mol_com_cluster_fit.xtc -s ${gmxdir}/md.tpr -n ${gmxdir}/index_test.ndx -b 100000 \
            -o ${outdir}/hbonds_gmx.ndx -pf yes -num ${outdir}/hbonds_gmx_num.xvg << EOL
Chain_A
Chain_B
EOL