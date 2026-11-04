"""Find and remove residues of a PDB file by residue name."""

from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdbio


def _is_atom_line(line: str) -> bool:
    return line.startswith(("ATOM  ", "HETATM"))


def _resname(line: str) -> str:
    return line[17:20].strip()


def remove_residues(pdb: Path, names, out: Path) -> list[str]:
    """
    Write a copy of a PDB file without the residues whose name is in ``names``.

    Every other line is copied unchanged. Returns the removed residues as
    "NAME chain resSeq", in file order.
    """
    names = set(names)
    removed, kept = [], []
    for line in Path(pdb).read_text().splitlines(keepends=True):
        if _is_atom_line(line) and _resname(line) in names:
            # One label per residue (name, chain, number and insertion code).
            label = f"{_resname(line)} {line[21].strip()}{line[22:27].strip()}"
            if label not in removed:
                removed.append(label)
            continue
        kept.append(line)
    Path(out).write_text("".join(kept))
    return removed


def polymer_residues(pdb: Path, names) -> list[str]:
    """
    Residues with a name in ``names`` that are amino acids or nucleotides.

    The classification uses the PDB Chemical Component Dictionary, so modified
    residues (e.g. MSE, SEP) count as amino acids. Returns "NAME chain resSeq" labels.
    """
    atoms = pdbio.PDBFile.read(str(pdb)).get_structure(model=1)
    polymer = struc.filter_amino_acids(atoms) | struc.filter_nucleotides(atoms)
    selected = atoms[polymer & np.isin(atoms.res_name, list(names))]
    labels = []
    for res in struc.residue_iter(selected):
        label = f"{res.res_name[0]} {res.chain_id[0]}{res.res_id[0]}"
        if label not in labels:
            labels.append(label)
    return labels
