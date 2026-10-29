#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Remove 5'-terminal phosphates from DNA/RNA strands before LEaP.

AMBER's 5'-terminal templates (G5, DA5, ...) describe a 5'-OH end, so a 5' phosphate
(P, OP1, OP2, OP3 on a strand's first nucleotide) is left untyped by LEaP and makes
saveAmberParm fail. This step writes <name>_trim5p.pdb: the input minus those atoms,
in every model, with everything else (headers, TER, formatting) unchanged. Strands that
start after an internal chain break are trimmed too, with a warning about the gap.
"""

import argparse
import sys
from pathlib import Path

import biotite.structure.io.pdb as pdbio

from scripts.utils.helpers import setup_logger
from scripts.utils.nucleic.infer import LINK_MAX
from scripts.utils.ter_phos import find_5prime_phosphates, remove_atom_lines
from scripts.utils.pdbter import ter_residues


def main(logfile: Path, pdb: Path, link_max: float = LINK_MAX) -> Path:
    """
    Trim 5'-terminal phosphates of every DNA/RNA strand.

    Arguments:
        logfile  : Path to the log file.
        pdb      : Input PDB file (not modified).
        link_max : Longest O3'-P distance (Å) still counted as a backbone link; a longer
                   one is a chain break, and the nucleotide after it starts a new strand
                   (default: 4.0).

    Returns:
        Path of the trimmed file, <name>_trim5p.pdb next to the input.
    """
    logger = setup_logger("prep", logfile)
    pdb = Path(pdb)
    out = pdb.with_name(f"{pdb.stem}_trim5p.pdb")

    ##### CHECK INPUT #####
    if not pdb.is_file():
        logger.error(f"PDB file {pdb} not found.")
        sys.exit(1)
    if link_max <= 0:
        logger.error("link_max must be positive.")
        sys.exit(1)

    ##### FIND 5'-TERMINAL PHOSPHATES (decided on model 1, applied to all models) #####
    pdb_file = pdbio.PDBFile.read(str(pdb))
    atoms = pdb_file.get_structure(model=1, altloc="occupancy")
    n_models = pdb_file.get_model_count()
    trims = find_5prime_phosphates(atoms, link_max, ter_after=ter_residues(pdb, model=1))

    ##### WRITE #####
    # newline="" keeps the original line endings, so untouched lines stay byte-identical.
    with open(pdb, newline="") as f:
        lines = f.readlines()
    new_lines, n_removed = remove_atom_lines(lines, trims)
    with open(out, "w", newline="") as f:
        f.writelines(new_lines)

    if not trims:
        logger.debug(f"No 5'-terminal phosphates found in {pdb.name}; {out.name} is an unchanged copy.")
        return out

    for t in trims:
        where = (f"after a chain break (O3'-P > {link_max} Å)" if t["after_break"]
                 else "first nucleotide of the chain")
        logger.info(f"Trimmed 5'-phosphate of {t['res_name']}{t['res_id']}{t['ins_code']} "
                    f"(chain {t['chain']}, {where}): removed {', '.join(t['atoms'])}")
        if t["after_break"]:
            logger.warning(f"Chain {t['chain']} has a gap before {t['res_name']}{t['res_id']}"
                           f"{t['ins_code']}; its dangling 5'-phosphate was removed, but the gap itself is not repaired. Review the structure.")
    return out
