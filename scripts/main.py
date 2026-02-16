import argparse
from scripts.pipeconfig import Pipeline

def main():
    parser = argparse.ArgumentParser(description="Pipeline for preparing and running MD simulations.")
    parser.add_argument('--params', dest='paramsfile', type=str, default='scripts/params.yml', help='Path to the YAML parameter file.')
    parser.add_argument('--prep', action='store_true', help='Run preparation step.')
    parser.add_argument('--md', action='store_true', help='Run MD simulation step.')
    parser.add_argument('--post', action='store_true', help='Run post-processing step.')
    parser.add_argument('--analysis', action='store_true', help='Run analysis step.')
    parser.add_argument('--rules', dest='rules', nargs='+', default=None, help='Snakefile rules to be used in analysis step.')
    parser.add_argument('--unlock', dest='unlock', action='store_true', help='Unlock snakemake workflow if it is locked.')
    args = parser.parse_args()

    # Build array of steps to run
    steps = []
    if args.prep:
        steps.append("prep")
    if args.md:
        steps.append("md")
    if args.post:
        steps.append("post")
    if args.analysis:
        steps.append("analysis")

    pipeline = Pipeline(paramsfile=args.paramsfile, 
                        steps=steps, 
                        rules=args.rules, 
                        unlock=args.unlock)
    pipeline.run_all()


if __name__ == "__main__":
    main()