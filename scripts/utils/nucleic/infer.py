"""Infer Watson-Crick pairs along duplex registers.

Pairs observed by geometry define *registers*: for two bases at positions i and j
(0-based index within their chain), the register is (chain1, chain2, i + j). Along
an antiparallel duplex or hairpin stem every pair has the same i + j.

For each register with at least ``min_observed`` observed pairs:

* gaps between observed pairs are filled with Watson-Crick pairs (see
  :func:`is_wc_pair`: A-T, G-C, and A-U only with RNA uracil); others are skipped;
* the register is extended outwards (towards the strand ends) while the bases
  are complementary; the first mismatch, used base, chain break or strand end stops it;
* hairpins (both bases on one chain) are never extended inwards: no pair may have
  a smaller span (j - i) than the innermost observed pair, so loops stay unpaired.

Inferred pairs never reuse a base that is already paired. "Next base" follows the
covalent O3'-P link (<= LINK_MAX, tolerant of stretched bonds in medium-resolution
models such as 1BF5 B1010-B1011 at 3.6 A), so a missing residue breaks the walk instead of shifting the
register.
"""

import numpy as np
import pandas as pd
import biotite.structure as struc

LINK_MAX = 4.0  # Angstrom; a bond is ~1.6 A, a missing nucleotide leaves a much larger gap


def is_wc_pair(letter1, rna1, letter2, rna2):
    """True if two bases form a canonical Watson-Crick pair that should be reported.

    Parameters
    ----------
    letter1, letter2 : str or None
        Parent base letters from ``map_nucleotide`` ("A", "C", "G", "T", "U"),
        or None if the base could not be identified.
    rna1, rna2 : bool
        True if the residue is a ribonucleotide (has an O2' atom).
    """
    if letter1 is None or letter2 is None:
        return False
    pair = frozenset((letter1, letter2))
    if pair in (frozenset("AT"), frozenset("GC")):
        return True
    if pair == frozenset("AU"):
        # Only RNA uracil counts; deoxyuridine in DNA is deliberately excluded.
        return rna1 if letter1 == "U" else rna2
    return False


def nucleic_residues(atoms, chains):
    """Residues of each nucleic-acid chain in file order, plus covalent links.

    Returns
    -------
    dict
        ``{chain: {"residues": [dict, ...], "links": [bool, ...]}}``. Each residue dict
        has ``chain``, ``res_id``, ``ins_code``, ``res_name``, ``letter`` (parent base
        from ``map_nucleotide``, or None) and ``is_rna`` (O2' present). ``links[k]`` is
        True when residue k's O3' is bonded to residue k+1's P.
    """
    nuc = atoms[struc.filter_nucleotides(atoms) & np.isin(atoms.chain_id, chains)]
    out = {}
    for chain_id in pd.unique(nuc.chain_id):
        chain = nuc[nuc.chain_id == chain_id]
        residues, o3, p = [], [], []
        for res in struc.residue_iter(chain):
            residues.append({"chain": chain_id, "res_id": int(res.res_id[0]),
                             "ins_code": res.ins_code[0], "res_name": res.res_name[0],
                             "letter": struc.map_nucleotide(res)[0],
                             "is_rna": bool((res.atom_name == "O2'").any())})
            o3.append(res.coord[res.atom_name == "O3'"])
            p.append(res.coord[res.atom_name == "P"])
        links = [len(o3[k]) > 0 and len(p[k + 1]) > 0
                 and np.linalg.norm(o3[k][0] - p[k + 1][0]) <= LINK_MAX
                 for k in range(len(residues) - 1)]
        out[chain_id] = {"residues": residues, "links": links}
    return out


def infer_pairs(letters, links, observed, min_observed=2, rna=None):
    """Infer pairs along the registers defined by the observed pairs.

    Parameters
    ----------
    letters : dict
        ``{chain: [letter or None, ...]}`` in chain order.
    links : dict
        ``{chain: [bool, ...]}``; ``links[c][k]`` is True if residue k links to k+1.
    observed : list
        Observed pairs as ``((chain1, i), (chain2, j))`` with 0-based indices.
    min_observed : int
        Observed pairs a register needs before it is used for inference.
    rna : dict, optional
        ``{chain: [bool, ...]}``, True for ribonucleotides. Default: all DNA.

    Returns
    -------
    inferred : list
        Inferred pairs in the same format, ordered so that ``(chain1, i) < (chain2, j)``.
    registers : list of dict
        ``chain1``, ``chain2``, ``sum``, ``n_observed``, ``n_inferred`` per register,
        most observed pairs first.
    """
    def ordered(a, b):
        return (a, b) if a <= b else (b, a)

    def is_rna(chain, k):
        return bool(rna and rna[chain][k])

    groups = {}
    for a, b in observed:
        (c1, i), (c2, j) = ordered(a, b)
        groups.setdefault((c1, c2, i + j), []).append(i)

    used = {base for pair in observed for base in pair}
    inferred, registers = [], []

    for (c1, c2, s), xs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        reg = {"chain1": c1, "chain2": c2, "sum": s, "n_observed": len(xs), "n_inferred": 0}
        registers.append(reg)
        if len(xs) < min_observed:
            continue
        xs = sorted(xs)
        hairpin = c1 == c2
        n1, n2 = len(letters[c1]), len(letters[c2])

        def linked(c, a, b):
            """True if every residue from index a to b is covalently linked to the next."""
            lo, hi = min(a, b), max(a, b)
            return all(links[c][k] for k in range(lo, hi))

        def try_add(i):
            j = s - i
            if not (0 <= i < n1 and 0 <= j < n2) or (hairpin and i >= j):
                return False
            if (c1, i) in used or (c2, j) in used:
                return False
            if not is_wc_pair(letters[c1][i], is_rna(c1, i), letters[c2][j], is_rna(c2, j)):
                return False
            used.update({(c1, i), (c2, j)})
            inferred.append(((c1, i), (c2, j)))
            reg["n_inferred"] += 1
            return True

        # 1. Fill gaps between observed pairs (skip mismatches, don't stop).
        for a, b in zip(xs, xs[1:]):
            if b - a > 1 and linked(c1, a, b) and linked(c2, s - b, s - a):
                for i in range(a + 1, b):
                    try_add(i)

        # 2. Extend outwards: i decreases, j increases (span grows for hairpins).
        i = xs[0]
        while i - 1 >= 0 and s - i + 1 < n2 and links[c1][i - 1] and links[c2][s - i]:
            i -= 1
            if not try_add(i):
                break

        # 3. Two strands only: also extend the other way (i increases, j decreases).
        #    For a hairpin this direction leads into the loop, so it is never taken.
        if not hairpin:
            i = xs[-1]
            while i + 1 < n1 and s - i - 1 >= 0 and links[c1][i] and links[c2][s - i - 1]:
                i += 1
                if not try_add(i):
                    break

    return [ordered(a, b) for a, b in inferred], registers
