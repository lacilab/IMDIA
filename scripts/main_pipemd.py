import argparse
import yaml
import prep01.prep
import md02.md

def main(datadir, filename, gmxdir, repn,
        with_dna = False, with_phosphate = False, range_a = None, range_b = None, ffs = None,
        ntmpi = None, ntomp = None, ngpu = None, gpulist = None, 
        em = None, nvt = None, npt = None, md = None):

        if args.prep:
            prep01.prep.main(datadir, filename, with_dna, with_phosphate, range_a, range_b, ffs)
        if args.md:
            md02.md.main(gmxdir, repn, with_dna, with_phosphate, ntmpi, ntomp, ngpu, gpulist, em, nvt, npt, md)


parser = argparse.ArgumentParser(description="Pipeline for preparing and running MD simulations.")
parser.add_argument('--prep', action='store_true', help='Run preparation step.')
parser.add_argument('--md', action='store_true', help='Run MD simulation step.')
args = parser.parse_args()

with open('main_params.yml') as f:
    params = yaml.safe_load(f)
    print("Loaded parameters from main_params.yml:", params)

main(**params)