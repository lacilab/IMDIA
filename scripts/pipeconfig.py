from dataclasses import dataclass
from pathlib import Path
import sys
import yaml


@dataclass
class Paths:
    datadir: Path
    mutdir: Path
    pdbfile: str
    gmxdir: Path
    paramspath: Path

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
    em: Path
    nvt: Path
    npt: Path
    md: Path
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
        paramspath = Path(self.paramsfile)
        self._validate_file(paramspath, "Parameter")
        mutdir = Path(pdbpath.parent)
        self._validate_dir(mutdir, "Mutation")
        datadir = Path(mutdir.parent.parent)
        self._validate_dir(datadir, "Data")
        gmxdir = mutdir / "gromacs_common" # Don't validate, will be created during prep step
        pdbfile = pdbpath.name 
        self._validate_file(mutdir / pdbfile, "In the mutation directory, PDB")
        return Paths(mutdir=mutdir, datadir=datadir, pdbfile=pdbfile,
                     gmxdir=gmxdir, paramspath=paramspath)
    
    def _validate_file(self, file, role) -> None:
        if not file.is_file():
            print(f"ERROR: {role} file '{file}' not found.")
            sys.exit(1)
        else:
            print(f"{role} file '{file}' found.")

    def _to_path(self, value):
        """Convert str/None to Path, or pass through existing Path."""
        if value is None:
            return None
        if isinstance(value, Path):
            return value
        return Path(value)

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

        pdbpath = self._to_path(params.get("pdbpath"))
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
            ntmpi=params.get("ntmpi", 0),
            ntomp=params.get("ntomp", 0),
            nb=params.get("nb", "auto"),
            pme=params.get("pme", "auto"),
            npme=params.get("npme", -1),
            ngpu=params.get("ngpu", None),
            gputasks=params.get("gputasks", None),
            em=self._to_path(params.get("em")),
            nvt=self._to_path(params.get("nvt")),
            npt=self._to_path(params.get("npt")),
            md=self._to_path(params.get("md")),
            snakecores=params.get("snakecores", None),
        )

        return pdbpath, system, run

    def setup(self) -> None:
        # Create log directory if it doesn't exist
        # logdir = self.paths.mutdir / "logs"
        # logdir.mkdir(exist_ok=True)
        # print(f"Log directory '{logdir}' is ready.")
        pass

    def prep(self) -> None:
        import scripts.prep01.prep_wrap as prep_wrap
        prep_wrap.main(self)

    def md(self) -> None:
        import scripts.md02.md_wrap as md_wrap
        md_wrap.main(self)

    def post(self) -> None:
        import scripts.post03.post_wrap as post_wrap
        post_wrap.main(self)

    def analysis(self) -> None:
        import scripts.analysis04.analysis_wrap as analysis_wrap
        analysis_wrap.main(self)
        
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
