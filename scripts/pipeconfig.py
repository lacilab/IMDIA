from dataclasses import dataclass
from pathlib import Path
import sys
import os
import yaml

import scripts.prep01.prep_wrap as prep_wrap
import scripts.md02.md_wrap as md_wrap
import scripts.post03.post_wrap as post_wrap
import scripts.analysis04.analysis_wrap as analysis_wrap

@dataclass
class Paths:
    datadir: Path
    mutdir: Path
    pdbfile: Path
    gmxdir: Path
    paramsfile: Path

@dataclass
class SystemConfig:
    repn: list[int]
    dna: bool
    phos: bool
    range_a: list[int]
    range_b: list[int]
    ffs: list[str]

@dataclass
class RunConfig:
    ntmpi: int
    ntomp: int
    nb: str
    pme: str
    npme: int
    ngpu: int
    gputasks: int
    em: bool
    nvt: bool
    npt: bool
    md: bool
    snakecores: int

class Pipeline:
    def __init__(self, paramsfile: Path | str, steps: list[str], rules: list[str] = None, unlock: bool = False):        
        self.steps = steps
        self.rules = rules
        self.unlock = unlock
        self.paramsfile = Path(paramsfile)
        self.pdbpath, self.system, self.run = self._get_params()
        self.paths = self._derive_paths()

    def _derive_paths(self) -> 'Paths':
        pdbpath = Path(self.pdbpath)
        self._validate_file(pdbpath, "PDB")
        paramsfile = self.paramsfile
        self._validate_file(paramsfile, "Parameter")
        mutdir = pdbpath.parent
        self._validate_dir(mutdir, "Mutation directory")
        datadir = mutdir.parent.parent
        self._validate_dir(datadir, "Data directory")
        gmxdir = mutdir / "gromacs_common" # Don't validate, will be created during prep step
        pdbfile = pdbpath.name 
        self._validate_file(mutdir / pdbfile, "PDB in mutation directory")
        return Paths(mutdir=mutdir, datadir=datadir, pdbfile=pdbfile,
                     gmxdir=gmxdir, paramsfile=paramsfile)
    
    def _validate_file(self, file, role) -> None:
        if not file.is_file():
            print(f"ERROR: {role} file '{file}' not found.")
            sys.exit(1)
        else:
            print(f"{role} file '{file}' found.")

    def _validate_dir(self, file, role) -> None:
        if not file.is_dir():
            print(f"ERROR: {role} directory '{file}' not found.")
            sys.exit(1)
        else:
            print(f"{role} directory '{file}' found.")
    
    
    def _get_params(self) -> tuple[str, SystemConfig, RunConfig]:
        if not self.paramsfile.is_file():
            print(f"ERROR: Parameter file '{self.paramsfile}' not found.")
            sys.exit(1)

        # Load parameters from YAML file
        with open(self.paramsfile) as f:
            params = yaml.safe_load(f)

        pdbpath = params.get("pdbpath", None)
        if pdbpath is None:
            print("ERROR: 'pdbpath' not specified in parameter file.")
            sys.exit(1)

        system = SystemConfig(
            repn=params.get("repn", []),
            dna=params.get("dna", False),
            phos=params.get("phos", False),
            range_a=params.get("range_a", [0, -1]),
            range_b=params.get("range_b", [0, -1]),
            ffs=params.get("ffs", None)
        )

        run = RunConfig(
            ntmpi=params.get("ntmpi", None),
            ntomp=params.get("ntomp", None),
            nb=params.get("nb", None),
            pme=params.get("pme", None),
            npme=params.get("npme", None),
            ngpu=params.get("ngpu", None),
            gputasks=params.get("gputasks", None),
            em=params.get("em", None),
            nvt=params.get("nvt", None),
            npt=params.get("npt", None),
            md=params.get("md", None),
            snakecores=params.get("snakecores", None),
        )

        return pdbpath, system, run

    def setup(self) -> None:
        # Create log directory if it doesn't exist
        logdir = self.paths.mutdir / "logs"
        os.makedirs(logdir, exist_ok=True)
        print(f"Log directory '{logdir}' is ready.")

    def prep(self) -> None:
        prep_wrap.main(self.paths.mutdir, 
                       self.paths.pdbfile, 
                       self.system.dna, 
                       self.system.phos, 
                       self.system.range_a, 
                       self.system.range_b, 
                       self.system.ffs)

    def md(self) -> None:
        for n in self.system.repn:
                md_wrap.main(self.paths.mutdir, 
                             self.paths.gmxdir, 
                             n, 
                             self.system.dna, 
                             self.system.phos, 
                             self.run.ntmpi, 
                             self.run.ntomp, 
                             self.run.nb, 
                             self.run.pme, 
                             self.run.npme, 
                             self.run.ngpu, 
                             self.run.gputasks, 
                             self.run.em, 
                             self.run.nvt, 
                             self.run.npt, 
                             self.run.md)

    def post(self) -> None:
        for n in self.system.repn:
            post_wrap.main(self.paths.mutdir, 
                           n, 
                           self.system.dna, 
                           self.system.phos)

    def analysis(self) -> None:
        if self.rules is None:
            analysis_wrap.main(self.paths.datadir, 
                               self.paths.paramsfile, 
                               self.run.snakecores, 
                               self.unlock)
        else:
            analysis_wrap.main(self.paths.datadir, 
                               self.paths.paramsfile, 
                               self.run.snakecores, 
                               self.rules, 
                               self.unlock)


    def run_all(self) -> None:
        self.setup()
        print(f"Running pipeline with steps: {self.steps}")
        if "prep" in self.steps:
            self.prep()
        if "md" in self.steps:
            self.md()
        if "post" in self.steps:
            self.post()
        if "analysis" in self.steps:
            self.analysis()
