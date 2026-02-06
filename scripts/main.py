import argparse
from datetime import datetime
import os
import sys
import yaml
import scripts.prep01.prep_wrap as prep_wrap
import scripts.md02.md_wrap as md_wrap
import scripts.post03.post_wrap as post_wrap
import scripts.analysis04.analysis_wrap as analysis_wrap
from scripts.utils.helpers import setup_logger


def pipeline(mutdir: str, datadir: str, pdbfile: str, gmxdir: str, args: argparse.Namespace, params: dict = None):

        # Defaults
        params = params or {}
        repn = params.get("repn", [])
        dna = params.get("dna", False) 
        phos = params.get("phos", False) 
        range_a = params.get("range_a", None)
        range_b = params.get("range_b", None)
        ffs = params.get("ffs", None)
        ntmpi = params.get("ntmpi", None)
        ntomp = params.get("ntomp", None)
        nb = params.get("nb", None) 
        pme = params.get("pme", None)
        npme = params.get("npme", None) 
        ngpu = params.get("ngpu", None)
        gputasks = params.get("gputasks", None) 
        em = params.get("em", None)
        nvt = params.get("nvt", None)
        npt = params.get("npt", None)
        md = params.get("md", None)
        snakecores = params.get("snakecores", None)

        if args.prep:
            prep_wrap.main(mutdir, pdbfile, dna, phos, range_a, range_b, ffs)
        if args.md:
            for n in repn:
                md_wrap.main(mutdir, gmxdir, n, dna, phos, ntmpi, ntomp, nb, pme, npme, ngpu, gputasks, em, nvt, npt, md)
        if args.post:
            for n in repn:
                post_wrap.main(mutdir, n, dna, phos)
        if args.analysis and args.rules is None:
            analysis_wrap.main(datadir, args.params, snakecores, ffs)
        elif args.analysis and args.rules is not None:
            analysis_wrap.main(datadir, args.params, snakecores, args.rules, ffs)
        
def get_dirs(pdbpath: str):
    if not os.path.isfile(pdbpath):
        print(f"ERROR: PDB file '{pdbpath}' not found.")
        sys.exit(1)
    mutdir = os.path.dirname(pdbpath)
    datadir = os.path.dirname(os.path.dirname(mutdir))
    gmxdir = os.path.join(mutdir, "gromacs_common")
    pdbfile = os.path.basename(pdbpath)
    return mutdir, datadir, pdbfile, gmxdir

def main():
    parser = argparse.ArgumentParser(description="Pipeline for preparing and running MD simulations.")
    parser.add_argument('--params', dest='params', type=str, default='scripts/params.yml', help='Path to the YAML parameter file.')
    parser.add_argument('--prep', action='store_true', help='Run preparation step.')
    parser.add_argument('--md', action='store_true', help='Run MD simulation step.')
    parser.add_argument('--post', action='store_true', help='Run post-processing step.')
    parser.add_argument('--analysis', action='store_true', help='Run analysis step.')
    parser.add_argument('--rules', dest='rules', nargs='+', default=None, help='Snakefile rules to be used in analysis step.')
    args = parser.parse_args()

    # Load parameters from YAML file
    if os.path.isfile(args.params):
        with open(args.params) as f:
            params = yaml.safe_load(f)
            print("Loaded parameters from main_params.yml:", params)
    else:
        print(f"ERROR: Parameter file '{args.params}' not found.")
        sys.exit(1)

    mutdir, datadir, pdbfile, gmxdir = get_dirs(params["pdbpath"])

    # Initialize log directory
    if not os.path.isdir(mutdir):
        print(f"ERROR: mutdir '{mutdir}' does not exist or is not a directory.")
        sys.exit(1)
    elif not os.path.isdir(mutdir + "/logs"):
        log_dir = os.path.join(mutdir, "logs")
        os.makedirs(log_dir, exist_ok=True)

    pipeline(mutdir, datadir, pdbfile, gmxdir, args, params)


if __name__ == "__main__":
    main()