"""Find and remove 5'-terminal phosphates of DNA/RNA strands.

AMBER's 5'-terminal residue templates (G5, DA5, ...) describe a 5'-OH end. A 5'
phosphate in the input (P, OP1, OP2, OP3 on a strand's first nucleotide) has no
template, so LEaP leaves those atoms untyped and saveAmberParm fails.

Strands are found from the covalent backbone, independent of residue names: a strand
starts at a chain's first nucleotide, after a non-nucleotide residue, or where the
O3'(previous)-P(this) distance exceeds ``link_max`` (a chain break). The output file is
the input with exactly the trimmed atoms' lines removed (all models), so headers, TER
records and formatting are preserved.
"""

import numpy as np
import biotite.structure as struc

from .naming import SUGAR_ATOMS
from .nucleic.infer import LINK_MAX
from .pdbter import residue_key

PHOSPHATE_ATOMS = {"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P"}
OH_BOND_MAX = 1.25   # Å; hydrogens this close to a removed oxygen go with it


def find_5prime_phosphates(atoms, link_max=LINK_MAX, ter_after=()):
    """Phosphate atoms to remove from the first nucleotide of every strand.

    Parameters
    ----------
    atoms : AtomArray
        One model (the decision applies to every model of the file).
    link_max : float
        Longest O3'-P distance (Angstrom) still counted as a backbone link.
    ter_after : set, optional
        ``(chain, res_id, ins_code)`` of residues followed by a TER record. The next residue starts a new strand.

    Returns
    -------
    list of dict
        One entry per trimmed residue: ``chain``, ``res_id``, ``ins_code``,
        ``res_name``, ``atoms`` (names to remove) and ``after_break`` (True if the strand
        starts after an internal chain break, rather than at a chain start or after TER).
    """

    trims = []
    for chain_id in dict.fromkeys(atoms.chain_id):
        chain = atoms[atoms.chain_id == chain_id]
        starts = struc.get_residue_starts(chain, add_exclusive_stop=True)
        prev = None                                   # previous nucleotide, if linked-eligible
        for k, (a, b) in enumerate(zip(starts[:-1], starts[1:])):
            res = chain[a:b]
            names = set(res.atom_name)
            if k > 0 and residue_key(chain, starts[k - 1]) in ter_after:
                prev = None                           # TER: strand ends here
            if not SUGAR_ATOMS <= names:              # not a nucleotide: next one starts a strand
                prev = None
                continue
            strand_start = prev is None or not _linked(prev, res, link_max)
            phosphate = names & PHOSPHATE_ATOMS
            if strand_start and phosphate:
                remove = phosphate | _hydrogens_on(res, phosphate)
                trims.append({"chain": chain_id, "res_id": int(res.res_id[0]),
                              "ins_code": res.ins_code[0], "res_name": res.res_name[0],
                              "atoms": sorted(remove), "after_break": prev is not None})
            prev = res
    return trims


def _linked(prev, res, link_max):
    o3 = prev.coord[prev.atom_name == "O3'"]
    p = res.coord[res.atom_name == "P"]
    return len(o3) > 0 and len(p) > 0 and np.linalg.norm(o3[0] - p[0]) <= link_max


def _hydrogens_on(res, phosphate):
    """Names of hydrogens whose nearest heavy atom is a removed phosphate atom.

    Judged by geometry, not by name: e.g. 2KOC labels its terminal-phosphate hydrogen
    "HO5'", but it sits 0.96 A from OP3 and 2.8 A from O5'.
    """
    is_h = np.isin(res.element, ["H", "D"])
    heavy, hydrogens = res[~is_h], res[is_h]
    if not hydrogens.array_length() or not heavy.array_length():
        return set()
    dist = np.linalg.norm(hydrogens.coord[:, None, :] - heavy.coord[None, :, :], axis=2)
    nearest = dist.argmin(axis=1)
    on_phosphate = (np.isin(heavy.atom_name[nearest], list(phosphate))
                    & (dist.min(axis=1) <= OH_BOND_MAX))
    return set(hydrogens.atom_name[on_phosphate])


def remove_atom_lines(lines, trims):
    """Drop the ATOM/HETATM lines of the trimmed atoms (every model) and fix CONECT records.

    Returns the new lines and the number of atom lines removed.
    """
    remove = {(t["chain"], t["res_id"], t["ins_code"], name) for t in trims for name in t["atoms"]}
    out, removed_serials, n_removed = [], set(), 0
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            key = (line[21], _int(line[22:26]), line[26].strip(), line[12:16].strip())
            if key in remove:
                n_removed += 1
                serial = _int(line[6:11])
                if serial is not None:
                    removed_serials.add(serial)
                continue
        out.append(line)
    if removed_serials:
        out = [_fix_conect(line, removed_serials) if line.startswith("CONECT") else line
               for line in out]
        out = [line for line in out if line is not None]
    return out, n_removed


def _fix_conect(line, removed):
    """CONECT record without references to removed serials; None if nothing is left."""
    fields = [line[k:k + 5] for k in range(6, 31, 5)]
    serials = [_int(f) for f in fields if f.strip()]
    if None in serials:             # non-numeric (e.g. hybrid-36): leave untouched
        return line
    if serials[0] in removed:
        return None
    partners = [s for s in serials[1:] if s not in removed]
    if not partners:
        return None
    ending = line[len(line.rstrip("\r\n")):] or "\n"     # keep the file's line ending
    return "CONECT" + "".join(f"{s:5d}" for s in [serials[0]] + partners) + ending


def _int(text):
    try:
        return int(text)
    except ValueError:
        return None
