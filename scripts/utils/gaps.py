"""Detect gaps in protein and nucleic-acid chains (with biotite).

Three kinds of finding, per chain, on polymer residues (amino acids and nucleotides;
water, ions, ligands and hydrogens are ignored):

* ``chain_break`` / ``backbone_distortion``: a backbone bond (C-N, O3'-P, or one inside
  a residue) outside [min_len, max_len], from ``biotite.structure.check_backbone_continuity``.
  Between residues it is a chain break (missing residues or a stretched link);
  inside a residue it points to a missing or misplaced backbone atom.
* ``missing_atoms``: heavy atoms of the residue's PDB Chemical Component Dictionary
  template (``biotite.structure.info.residue``) that are absent. Leaving atoms (OXT,
  OP3) and the phosphate of a strand's 5'-terminal nucleotide are not required.
* ``numbering_jump``: residue numbers that skip (``check_res_id_continuity``) while the
  backbone is continuous. Jumps that coincide with a chain break are merged into it and
  give its estimated size instead.
"""

import numpy as np
import pandas as pd
import biotite.structure as struc
from biotite.structure.info import residue as ccd_residue

from .naming import translate_residue_names
from .pdbter import residue_key

GAP_COLUMNS = ["type", "chain", "residue_before", "residue_after", "atoms", "distance_A",
               "est_missing_residues", "details"]
# Canonical backbone order, so the continuity check does not depend on the file's atom order.
BACKBONE_ORDER = {name: k for k, name in enumerate(["N", "CA", "C", "P", "O5'", "C5'", "C4'",
                                                    "C3'", "O3'"])}
LEAVING_ATOMS = {"OXT", "OP3"}
FIVE_PRIME_PHOSPHATE = {"P", "OP1", "OP2", "OP3"}
# Older atom names still written by some tools; applied to both the file and the template
# (C5M is the old thymine methyl name, but the correct name in e.g. 5MU).
ATOM_ALIASES = {"O1P": "OP1", "O2P": "OP2", "O3P": "OP3", "C5M": "C7"}


def polymer_atoms(atoms):
    """Heavy atoms of amino-acid and nucleotide residues, with PDB residue names."""
    atoms = translate_residue_names(atoms)
    atoms = atoms[~np.isin(atoms.element, ["H", "D"])]
    return atoms[struc.filter_amino_acids(atoms) | struc.filter_nucleotides(atoms)]


def _canonical(names):
    return {ATOM_ALIASES.get(n, n) for n in names}


def _label(chain, k):
    return f"{chain.res_name[k]}{chain.res_id[k]}{chain.ins_code[k]}"


def _canonical_order(chain):
    """Reorder atoms inside each residue so backbone atoms follow BACKBONE_ORDER."""
    residue_index = np.cumsum(np.isin(np.arange(chain.array_length()),
                                      struc.get_residue_starts(chain))) - 1
    rank = np.array([BACKBONE_ORDER.get(n, len(BACKBONE_ORDER)) for n in chain.atom_name])
    return chain[np.lexsort((rank, residue_index))]


def backbone_gaps(chain, min_len=1.2, max_len=1.8):
    """Chain breaks and backbone distortions in one chain."""
    chain = _canonical_order(chain)
    after = struc.check_backbone_continuity(chain, min_len, max_len)
    backbone = np.where(struc.filter_peptide_backbone(chain) | struc.filter_phosphate_backbone(chain))[0]
    rows = []
    for k in after:
        pos = np.searchsorted(backbone, k)
        if pos == 0:
            continue
        prev = backbone[pos - 1]
        dist = float(np.linalg.norm(chain.coord[k] - chain.coord[prev]))
        same = (chain.res_id[k], chain.ins_code[k]) == (chain.res_id[prev], chain.ins_code[prev])
        missing = int(chain.res_id[k] - chain.res_id[prev] - 1)
        rows.append({
            "type": "backbone_distortion" if same else "chain_break",
            "chain": chain.chain_id[0],
            "residue_before": _label(chain, prev), "residue_after": _label(chain, k),
            "atoms": f"{chain.atom_name[prev]}-{chain.atom_name[k]}",
            "distance_A": round(dist, 2),
            "est_missing_residues": None if same else max(missing, 0),
            "details": ("bond too short (atoms overlap)" if dist < min_len else
                        "bond too long" if same else
                        f"{max(missing, 0)} residue(s) missing by numbering" if missing > 0 else
                        "stretched link, numbering continuous"),
        })
    return rows


def missing_atoms(chain):
    """Residues with heavy atoms missing vs their CCD template, and residues not checked."""
    rows, not_checked = [], []
    starts = struc.get_residue_starts(chain, add_exclusive_stop=True)
    first_nucleotide = True
    for a, b in zip(starts[:-1], starts[1:]):
        res = chain[a:b]
        name = res.res_name[0]
        try:
            template = ccd_residue(name)
        except Exception:
            not_checked.append(_label(chain, a))
            continue
        expected = _canonical(template.atom_name[template.element != "H"]) - LEAVING_ATOMS
        is_nucleotide = bool(struc.filter_nucleotides(res).any())
        if is_nucleotide and first_nucleotide:
            expected -= FIVE_PRIME_PHOSPHATE
        first_nucleotide = first_nucleotide and not is_nucleotide
        absent = sorted(expected - _canonical(res.atom_name))
        if absent:
            rows.append({"type": "missing_atoms", "chain": chain.chain_id[0],
                         "residue_before": _label(chain, a), "residue_after": None,
                         "atoms": " ".join(absent), "distance_A": None,
                         "est_missing_residues": None,
                         "details": f"{len(absent)} of {len(expected)} heavy atoms missing"})
    return rows, not_checked


def numbering_jumps(chain, breaks):
    """Numbering skips not explained by a chain break."""
    broken = {(r["residue_before"], r["residue_after"]) for r in breaks if r["type"] == "chain_break"}
    starts = struc.get_residue_starts(chain)
    rows = []
    for k in struc.check_res_id_continuity(chain):
        prev = starts[np.searchsorted(starts, k) - 1]
        before, after = _label(chain, prev), _label(chain, k)
        if (before, after) in broken:
            continue
        rows.append({"type": "numbering_jump", "chain": chain.chain_id[0],
                     "residue_before": before, "residue_after": after, "atoms": None,
                     "distance_A": None, "est_missing_residues": None,
                     "details": f"numbering {chain.res_id[prev]} -> {chain.res_id[k]}, "
                                "backbone continuous (renumbering?)"})
    return rows


def chain_kind(chain):
    """'protein', 'DNA', 'RNA' or a mix, from the residues present."""
    starts = struc.get_residue_starts(chain, add_exclusive_stop=True)
    kinds = set()
    for a, b in zip(starts[:-1], starts[1:]):
        res = chain[a:b]
        if struc.filter_amino_acids(res).any():
            kinds.add("protein")
        else:  # nucleotide: RNA if the sugar has O2' (also for modified residues)
            kinds.add("RNA" if (res.atom_name == "O2'").any() else "DNA")
    return "/".join(sorted(kinds))


def segments(chain, ter_after):
    """Split one chain at TER records (a TER ends a strand, like a new chain ID)."""
    starts = struc.get_residue_starts(chain)
    cuts = [starts[k] for k in range(1, len(starts))
            if residue_key(chain, starts[k - 1]) in ter_after]
    bounds = [0] + cuts + [chain.array_length()]
    return [chain[a:b] for a, b in zip(bounds[:-1], bounds[1:])]


def diagnose(atoms, min_len=1.2, max_len=1.8, ter_after=()):
    """Find gaps in all polymer chains.

    ``ter_after`` holds residues followed by a TER record (see
    :func:`scripts.utils.pdbter.ter_residues`). A chain split by TER is analysed as
    separate segments, labelled "A/1", "A/2", ..., so the end of one strand and the start
    of the next are not reported as a gap.

    Returns
    -------
    gaps : pandas.DataFrame
        One row per finding, columns ``GAP_COLUMNS``.
    chains : pandas.DataFrame
        Per chain: kind, residues, chain breaks, backbone distortions, residues with
        missing atoms, numbering jumps, residues not checked.
    """
    polymer = polymer_atoms(atoms)
    rows, summary = [], []
    chains = []
    for chain_id in pd.unique(polymer.chain_id):
        parts = segments(polymer[polymer.chain_id == chain_id], ter_after)
        for n, part in enumerate(parts, start=1):
            if len(parts) > 1:
                part = part.copy()
                part.chain_id[:] = f"{chain_id}/{n}"
            chains.append(part)
    for chain in chains:
        chain_id = chain.chain_id[0]
        breaks = backbone_gaps(chain, min_len, max_len)
        missing, not_checked = missing_atoms(chain)
        jumps = numbering_jumps(chain, breaks)
        rows += breaks + missing + jumps
        summary.append({
            "chain": chain_id, "kind": chain_kind(chain),
            "residues": len(struc.get_residue_starts(chain)),
            "chain_breaks": sum(r["type"] == "chain_break" for r in breaks),
            "backbone_distortions": sum(r["type"] == "backbone_distortion" for r in breaks),
            "residues_missing_atoms": len(missing),
            "numbering_jumps": len(jumps),
            "not_checked": " ".join(not_checked),
        })
    return pd.DataFrame(rows, columns=GAP_COLUMNS), pd.DataFrame(summary)
