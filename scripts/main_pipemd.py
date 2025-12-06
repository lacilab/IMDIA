import argparse
from datetime import datetime
import os
import sys
import yaml
import scripts.prep01.prep_wrap as prep_wrap
import scripts.md02.md_wrap as md_wrap
import scripts.post03.post_wrap as post_wrap
from scripts.utils.helpers import setup_logger


def pipeline(args: argparse.Namespace, datadir: str, filename: str, gmxdir: str, repn: int,
        dna: bool = False, phos: bool = False, range_a: str = None, range_b: str = None, ffs: str = None,
        ntmpi: int = None, ntomp: int = None, nb: str = None, pme: str = None, npme: int = None, ngpu: int = None, gpulist: str = None, 
        em: str = None, nvt: str = None, npt: str = None, md: str = None):

        if args.prep:
            prep_wrap.main(datadir, filename, dna, phos, range_a, range_b, ffs)
        if args.md:
            md_wrap.main(datadir, gmxdir, repn, dna, phos, ntmpi, ntomp, nb, pme, npme, ngpu, gpulist, em, nvt, npt, md)
        if args.post:
            post_wrap.main(datadir, gmxdir, dna, phos)


def main():
    parser = argparse.ArgumentParser(description="Pipeline for preparing and running MD simulations.")
    parser.add_argument('--params', dest='params', type=str, default='scripts/params.yml', help='Path to the YAML parameter file.')
    parser.add_argument('--prep', action='store_true', help='Run preparation step.')
    parser.add_argument('--md', action='store_true', help='Run MD simulation step.')
    parser.add_argument('--post', action='store_true', help='Run post-processing step.')
    args = parser.parse_args()

    # Load parameters from YAML file
    if os.path.isfile(args.params):
        with open(args.params) as f:
            params = yaml.safe_load(f)
            print("Loaded parameters from main_params.yml:", params)
    else:
        print(f"ERROR: Parameter file '{args.params}' not found.")
        sys.exit(1)

    # Initialize log directory
    if not os.path.isdir(params["datadir"]):
        print(f"ERROR: datadir '{params['datadir']}' does not exist or is not a directory.")
        sys.exit(1)
    elif not os.path.isdir(params["datadir"] + "/logs"):
        log_dir = os.path.join(params["datadir"], "logs")
        os.makedirs(log_dir, exist_ok=True)
        
    pipeline(args, **params)


if __name__ == "__main__":
    main()