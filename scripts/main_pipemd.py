import argparse
import os
import sys
import yaml
import scripts.prep01.prep_wrap as prep_wrap
import scripts.md02.md_wrap as md_wrap
# from scripts.utils.helpers import error

def pipeline(args, datadir, filename, gmxdir, repn,
        dna = False, phos = False, range_a = None, range_b = None, ffs = None,
        ntmpi = None, ntomp = None, ngpu = None, gpulist = None, 
        em = None, nvt = None, npt = None, md = None):

        if args.prep:
            prep_wrap.main(datadir, filename, dna, phos, range_a, range_b, ffs)
        if args.md:
            md_wrap.main(datadir, gmxdir, repn, dna, phos, ntmpi, ntomp, ngpu, gpulist, em, nvt, npt, md)

def main():
    parser = argparse.ArgumentParser(description="Pipeline for preparing and running MD simulations.")
    parser.add_argument('--prep', action='store_true', help='Run preparation step.')
    parser.add_argument('--md', action='store_true', help='Run MD simulation step.')
    args = parser.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    yaml_path = os.path.join(here, "main_params.yml")
    with open(yaml_path) as f:
        params = yaml.safe_load(f)
        print("Loaded parameters from main_params.yml:", params)

    if not os.path.isdir(params["datadir"]):
        print(f"Error: datadir '{params['datadir']}' does not exist or is not a directory.")
        sys.exit(1)
    elif not os.path.isdir(params["datadir"] + "/logs"):
        log_dir = os.path.join(params["datadir"], "logs")
        os.makedirs(log_dir, exist_ok=True)

    pipeline(args, **params)

if __name__ == "__main__":
    main()