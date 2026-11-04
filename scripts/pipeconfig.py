from dataclasses import dataclass, replace
from datetime import datetime
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
    postpdb: Path | None    # post-simulation PDB used by the analysis stage (None without replicates)

@dataclass
class SystemConfig:
    repn: list[int]
    dna: bool | None    # None = detect from the structure
    phos: bool | None   # None = detect from the structure
    ffs: list[str]
    rna: bool | None = None  # None = detect from the structure
    protein: bool | None = None  # None = detect from the structure
    stop_on_gaps: bool = True    # stop before preparation on chain breaks or backbone distortions
    protonation: dict | None = None

@dataclass
class RunConfig:
    md_time: int
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
        # Pipeline-level log file for messages produced while setting up the run.
        self.logger = self._pipeline_logger()
        # Composition detected from the input structure; fills dna/rna/phos left on auto.
        self.composition = self._composition()
        self._resolve_flags()

    def _derive_paths(self) -> 'Paths':
        pdbpath = Path(self.pdbpath)
        self._validate_file(pdbpath, "PDB")
        paramspath = Path(self.paramsfile)
        self._validate_file(paramspath, "Parameter")
        mutdir = Path(pdbpath.parent)
        self._validate_dir(mutdir, "Mutation")
        datadir = Path(mutdir.parent.parent)
        self._validate_dir(datadir, "Data")
        gmxdir = mutdir / "gromacs_common" # Don't validate here, will be created during prep step
        # Post-simulation PDB of the first replicate (None without replicates).
        # Don't validate here, will be created during md step
        reps = self.system.repn
        postpdb = mutdir / f"rep{reps[0]}" / "gromacs" / "md_prot.pdb" if reps else None
        pdbfile = pdbpath.name 
        self._validate_file(mutdir / pdbfile, "In the mutation directory, PDB")
        return Paths(mutdir=mutdir, datadir=datadir, pdbfile=pdbfile,
                     gmxdir=gmxdir, paramspath=paramspath, postpdb=postpdb)
    
    def _validate_file(self, file, role) -> None:
        if not file.is_file():
            print(f"ERROR: {role} file '{file}' not found.")
            sys.exit(1)
        #else:
        #    print(f"{role} file '{file}' found.")

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
        #else:
        #    print(f"{role} directory '{file}' found.")
    
    
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
            dna=self._flag(params, "dna"),
            phos=self._flag(params, "phos"),
            rna=self._flag(params, "rna"),
            protein=self._flag(params, "protein"),
            ffs=params.get("ffs", None),
            stop_on_gaps=self._bool(params, "stop_on_gaps", True),
            protonation=params.get("protonation", None),
        )

        run = RunConfig(
            md_time=params.get("md_time", 10000), # Default duration 10 ns
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

    @staticmethod
    def _flag(params: dict, key: str) -> bool | None:
        """Read an optional true/false/auto key; None (auto) when missing or 'auto'."""
        value = params.get(key)
        if value is None or (isinstance(value, str) and value.lower() == "auto"):
            return None
        if isinstance(value, bool):
            return value
        print(f"ERROR: '{key}' in the parameter file must be true, false or auto (got {value!r}).")
        sys.exit(1)

    @staticmethod
    def _bool(params: dict, key: str, default: bool) -> bool:
        """Read an optional true/false key; the default when missing."""
        value = params.get(key)
        if value is None:
            return default
        if isinstance(value, bool):
            return value
        print(f"ERROR: '{key}' in the parameter file must be true or false (got {value!r}).")
        sys.exit(1)

    def _pipeline_logger(self):
        """Logger writing to mutdir/logs/pipeline_<timestamp>.log."""
        from scripts.utils.helpers import setup_logger

        logdir = self.paths.mutdir / "logs"
        logdir.mkdir(parents=True, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return setup_logger("pipeline", logdir / f"pipeline_{timestamp}.log")

    def _composition(self):
        """Detect the composition of the input PDB and keep a JSON record in mutdir.

        If a record from an earlier run exists, differences are reported as a warning
        (the input changed between runs). The record is then replaced by the new one.
        """
        from scripts.utils import composition as comp_mod

        record = self.paths.mutdir / comp_mod.RECORD_NAME
        comp = comp_mod.detect(self.paths.mutdir / self.paths.pdbfile)

        # Compare with the record of a previous run, if any.
        if record.is_file():
            try:
                previous = comp_mod.load(record)
                changes = comp_mod.differences(previous, comp)
            except Exception as err:
                previous = None
                changes = [f"previous record could not be read ({err})"]
            if changes:
                self.logger.warning(f"Composition differs from the previous record {record}:\n"
                                    + "\n".join(f"  - {change}" for change in changes))
            # A chain mapping stays valid only while the input file is unchanged.
            if previous is not None and previous.chain_mapping:
                if previous.source_sha256 == comp.source_sha256:
                    comp = replace(comp, chain_mapping=previous.chain_mapping)
                else:
                    self.logger.warning("The chain mapping in the previous record refers to "
                                        "a different input file and was removed; it is "
                                        "recreated by the analysis stage.")

        comp_mod.save(comp, record)
        self.logger.info(f"Detected composition: {comp.summary()} (record: {record})")
        return comp

    def _resolve_flags(self) -> None:
        """Set protein/dna/rna/phos from detection unless the parameter file sets them explicitly."""
        detected = {"protein": self.composition.has_protein,
                    "dna": self.composition.has_dna,
                    "rna": self.composition.has_rna,
                    "phos": self.composition.has_phospho}
        for key, found in detected.items():
            given = getattr(self.system, key)
            if given is None:
                # Not set by the user: use the detected value.
                setattr(self.system, key, found)
                self.logger.debug(f"{key}: {found} (detected)")
            elif given != found:
                # Explicit value that disagrees with the structure: keep it, but warn.
                self.logger.warning(
                    f"{key} is set to {given} in {self.paramsfile.name}, but the structure "
                    f"analysis found {found}. Using {given}. This is likely to cause errors "
                    "during the simulation if the forcefields required for the detected "
                    "molecule types are not included.")
            else:
                # Explicit value that agrees with the structure: confirmation only.
                self.logger.debug(f"{key}: {given} (set in {self.paramsfile.name}, "
                                  "matches detection)")

    def setup(self) -> None:
        # Create log directory if it doesn't exist
        # logdir = self.paths.mutdir / "logs"
        # logdir.mkdir(exist_ok=True)
        # print(f"Log directory '{logdir}' is ready.")
        pass

    def diag(self) -> None:
        import scripts.diag00.diag_wrap as diag_wrap
        diag_wrap.main(self)

    def prep(self) -> None:
        import scripts.prep01.prep_wrap as prep_wrap
        prep_wrap.main(self)

    def md(self, steps: list) -> None:
        import scripts.md02.md_wrap as md_wrap
        md_wrap.main(self, steps=steps)

    def post(self) -> None:
        import scripts.post03.post_wrap as post_wrap
        post_wrap.main(self)

    def analysis(self) -> None:
        # The analysis stage needs the post-simulation PDB.
        if self.paths.postpdb is None:
            print("ERROR: no replicates in 'repn' in the parameter file; the post-simulation "
                  "PDB (rep<n>/gromacs/md_prot.pdb) is needed for analysis.")
            sys.exit(1)
        self._validate_file(self.paths.postpdb, "Post-simulation PDB")
        import scripts.analysis04.analysis_wrap as analysis_wrap
        analysis_wrap.main(self)
        
    def run_all(self) -> None:
        self.setup()
        print(f"Running pipeline with steps: {self.steps}")
        # Diagnosis runs on its own (--diag) and always before preparation, once.
        if "diag" in self.steps or "prep" in self.steps:
            self.diag()
        if "prep" in self.steps:
            self.prep()
        # For items in self.steps that start with "md"
        md_steps = [step for step in self.steps if step.startswith("md")]
        if any(md_steps):
            self.md(md_steps)
        if "post" in self.steps:
            self.post()
        if "analysis" in self.steps:
            self.analysis()
