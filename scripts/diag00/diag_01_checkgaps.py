#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnose gaps in the protein and nucleic-acid chains of a PDB file.

Reports chain breaks (backbone bonds outside the expected length), residues with
missing heavy atoms, and residue-numbering jumps, then warns the user to review them.
Works with RCSB files and with modelled/processed files (AMBER/CHARMM residue names are
translated first). Missing residues at chain ends are not reported: they leave no
geometric trace.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

from scripts.utils.helpers import setup_logger
from scripts.utils.nucleic.io import load_structure
from scripts.utils.gaps import diagnose
from scripts.utils.pdbter import ter_residues

BLOCKING = ("chain_break", "backbone_distortion", "missing_atoms")
# Finding types that stop the pipeline when stop_on_gaps is true (the default).
STOPPING = ("chain_break", "backbone_distortion")

# Explanation and suggested fix per finding type, printed (only for the types found)
# between the gap warning and the list of findings.
ADVICE_HEADER = "What the findings mean and how to fix them (ChimeraX recommended):"
ADVICE = {
    "chain_break": """\
  [chain_break] Missing residues or a broken backbone link creating gaps. SEVERE. The preparation step might 
    not crash, but it may distort the structure near the gap.
    Fix: model the missing residues.
      Protein: ChimeraX Tools > Structure Editing > Model Loops (MODELLER; needs a
      free licence key).
      DNA/RNA: build the missing nucleotides in a modelling tool, or, if the break is
      real (two strands), end the strand with a TER record so it is treated as two.""",
    "backbone_distortion": """\
  [backbone_distortion] A backbone bond inside a residue is too long or too short,
    usually a misplaced or overlapping atom. MODERATE: it can cause very
    large forces in minimisation. Fix: inspect the residue in ChimeraX. Running Tools > Structure Editing > Dock Prep 
    often fixes it by picking a library rotamer that avoids clashes.""",
    "missing_atoms": """\
  [missing_atoms] Heavy atoms absent from a residue, usually a side chain not resolved
    in the experiment. MILD for surface residues, MODERATE in interfaces, binding sites
    or near a mutation. During preparation, LEaP will rebuild the atoms from its template geometry, without
    considering neighbouring atoms, so they may clash or point the wrong way.
    Fix: ChimeraX Tools > Structure Editing > Dock Prep,
    which picks a library rotamer that avoids clashes.""",
}
ADVICE_FOOTER = "  Re-run the diagnosis on the repaired structure before preparation."


def main(logfile: Path, pdb: Path, outdir: Path = None, model: int = 1,
         min_bond: float = 1.2, max_bond: float = 1.8, stop_on_gaps: bool = True):
    """
    Diagnose gaps in a PDB file and warn the user to review them.

    Arguments:
        logfile      : Path to the log file.
        pdb          : Input PDB file.
        outdir       : Folder for gaps.csv (default: the PDB's folder).
        model        : MODEL number to analyse (default: 1).
        min_bond     : Shortest accepted backbone bond (Å) (default: 1.2).
        max_bond     : Longest accepted backbone bond (Å) (default: 1.8).
        stop_on_gaps : Stop the pipeline (exit code 1) if chain breaks or backbone
                       distortions are found (default: True). Missing atoms never stop it.

    Returns:
        dict with the number of findings per type and the path of gaps.csv.
    """
    logger = setup_logger("diag", logfile)
    pdb = Path(pdb)
    outdir = Path(outdir) if outdir else pdb.parent

    ##### CHECK INPUT #####
    if not pdb.is_file():
        logger.error(f"PDB file {pdb} not found.")
        sys.exit(1)
    if not 0 < min_bond < max_bond:
        logger.error("Invalid bond window: ensure that 0 < min_bond < max_bond.")
        sys.exit(1)
    outdir.mkdir(parents=True, exist_ok=True)

    ##### DIAGNOSE #####
    try:
        atoms, info = load_structure(pdb, model=model)
    except ValueError as err:
        logger.error(str(err))
        sys.exit(1)
    for note in info["notes"]:
        logger.info(note)

    gaps, chains = diagnose(atoms, min_bond, max_bond, ter_after=ter_residues(pdb, model))
    if chains.empty:
        logger.warning(f"No protein or nucleic-acid residues found in {pdb.name}.")
        return {"gaps_csv": None}

    logger.info(f"Gap diagnosis of {pdb.name} (model {info['model_used']}):")
    for c in chains.itertuples():
        logger.debug(f"  chain {c.chain} ({c.kind}, {c.residues} residues): "
                    f"{c.chain_breaks} chain break(s), {c.backbone_distortions} backbone "
                    f"distortion(s), {c.residues_missing_atoms} residue(s) with missing atoms, "
                    f"{c.numbering_jumps} numbering jump(s)")
        if c.not_checked:
            logger.debug(f"    not checked (no residue template): {c.not_checked}")

    gaps_csv = outdir / f"{pdb.stem}_gaps.csv"
    gaps.to_csv(gaps_csv, index=False)
    counts = gaps["type"].value_counts().to_dict()
    blocking = gaps[gaps["type"].isin(BLOCKING)]
    jumps = gaps[gaps["type"] == "numbering_jump"]

    if blocking.empty:
        logger.info("No gaps found: the backbone is continuous and no heavy atoms are missing.")
    else:
        lines = [f"  [{g.type}] chain {g.chain}: {g.residue_before}"
                 + (f" -> {g.residue_after}" if pd.notna(g.residue_after) else "")
                 + (f" ({g.atoms})" if pd.notna(g.atoms) else "") + f": {g.details}"
                 for g in blocking.itertuples()]
        # Explanations for the finding types present, in BLOCKING order.
        found = set(blocking["type"])
        advice = ([ADVICE_HEADER] + [ADVICE[t] for t in BLOCKING if t in found]
                  + [ADVICE_FOOTER])
        logger.warning(f"{len(blocking)} gap(s) found in {pdb.name}.\n"
                       + "\n".join(advice) + "\nFindings:\n" + "\n".join(lines))
    if not jumps.empty:
        logger.info(f"{len(jumps)} residue-numbering jump(s) with a continuous backbone "
                    "(informational, often just renumbering): "
                    + "; ".join(f"chain {g.chain} {g.residue_before} -> {g.residue_after}"
                                for g in jumps.itertuples()))
    logger.info(f"Details written to {gaps_csv}")

    # Chain breaks and backbone distortions stop the pipeline unless stop_on_gaps is false.
    stopping = gaps[gaps["type"].isin(STOPPING)]
    if not stopping.empty and stop_on_gaps:
        logger.error(f"The pipeline stopped because chain breaks or backbone distortions were "
                     f"found ({len(stopping)}); stopping is the default. Repair the structure "
                     "(see the findings above) or set \"stop_on_gaps: false\" in the parameter "
                     "file to continue anyway; the gaps are then only reported.")
        sys.exit(1)
    if not stopping.empty:
        logger.warning(f"Continuing despite chain breaks or backbone distortions ({len(stopping)}) "
                       "because stop_on_gaps is false in the parameter file. Remove that setting "
                       "(default: true) to stop the pipeline when such gaps and distortions are found.")
    elif not blocking.empty:
        # Only missing atoms were found.
        logger.info("Continuing: missing atoms do not stop the pipeline (LEaP rebuilds them "
                    "during preparation); only chain breaks and backbone distortions do.")
    return {"gaps_csv": gaps_csv, **counts}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--logfile", type=Path, required=True, help="Log file")
    parser.add_argument("--pdb", type=Path, required=True, help="Input PDB file")
    parser.add_argument("--outdir", type=Path, default=None, help="Folder for the gaps CSV")
    parser.add_argument("--model", type=int, default=1, help="MODEL number (default 1)")
    parser.add_argument("--min-bond", dest="min_bond", type=float, default=1.2,
                        help="Shortest accepted backbone bond in Å (default 1.2)")
    parser.add_argument("--max-bond", dest="max_bond", type=float, default=1.8,
                        help="Longest accepted backbone bond in Å (default 1.8)")
    main(**vars(parser.parse_args()))
