import argparse
import os
import yaml
import scripts.prep01.prep_wrap as prep_wrap
import scripts.md02.md_wrap as md_wrap

def main(datadir, filename, gmxdir, repn,
        dna = False, phos = False, range_a = None, range_b = None, ffs = None,
        ntmpi = None, ntomp = None, ngpu = None, gpulist = None, 
        em = None, nvt = None, npt = None, md = None):

        if args.prep:
            prep_wrap.main(datadir, filename, dna, phos, range_a, range_b, ffs)
        if args.md:
            md_wrap.main(gmxdir, repn, dna, phos, ntmpi, ntomp, ngpu, gpulist, em, nvt, npt, md)

parser = argparse.ArgumentParser(description="Pipeline for preparing and running MD simulations.")
parser.add_argument('--prep', action='store_true', help='Run preparation step.')
parser.add_argument('--md', action='store_true', help='Run MD simulation step.')
args = parser.parse_args()

here = os.path.dirname(os.path.abspath(__file__))
yaml_path = os.path.join(here, "main_params.yml")
with open(yaml_path) as f:
    params = yaml.safe_load(f)
    print("Loaded parameters from main_params.yml:", params)

main(**params)