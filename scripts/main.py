import argparse
import os
from scripts.pipeconfig import Pipeline
from scripts.utils.helpers import VERBOSE_ENV

def main():
    parser = argparse.ArgumentParser(description="Pipeline for preparing and running MD simulations.")
    parser.add_argument('--params', dest='paramsfile', type=str, default='scripts/params.yml', help='Path to the YAML parameter file.')
    parser.add_argument('--prep', action='store_true', help='Run preparation step.')
    parser.add_argument('--md', action='store_true', help='Run full MD simulation (energy minimization, NVT and NPT equilibration, MD Production).')
    parser.add_argument('--md_em', action='store_true', help='Run energy minimization step.')
    parser.add_argument('--md_eq', action='store_true', help='Run equilibration step (NVT and NPT).')
    parser.add_argument('--md_prod', action='store_true', help='Run MD production step.')
    parser.add_argument('--post', action='store_true', help='Run post-processing step.')
    parser.add_argument('--analysis', action='store_true', help='Run analysis step.')
    parser.add_argument('--rules', dest='rules', nargs='+', default=None, help='Snakefile rules to be used in analysis step.')
    parser.add_argument('--unlock', dest='unlock', action='store_true', help='Unlock snakemake workflow if it is locked.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Also print DEBUG messages on the console; the log files always contain them.')
    args = parser.parse_args()

    # Loggers read this when they are first set up, so it is set before the pipeline starts.
    # Child processes (e.g. Snakemake rules) inherit it.
    if args.verbose:
        os.environ[VERBOSE_ENV] = "1"

    # Build array of steps to run
    steps_flags = [
    ("prep", args.prep),
    ("md", args.md),
    ("md_em", args.md_em),
    ("md_eq", args.md_eq),
    ("md_prod", args.md_prod),
    ("post", args.post),
    ("analysis", args.analysis),
    ]

    steps = [name for name, flag in steps_flags if flag]


    pipeline = Pipeline(paramsfile=args.paramsfile, 
                        steps=steps, 
                        rules=args.rules, 
                        unlock=args.unlock)
    pipeline.run_all()


if __name__ == "__main__":
    main()