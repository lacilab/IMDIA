"""Part 2: detect Watson-Crick base pairs (A-T, G-C; A-U with RNA uracil) with biotite.

Pipeline
--------
1. ``biotite.structure.base_pairs`` finds all base pairs (DSSR criteria,
   canonical and non-canonical).
2. ``map_nucleotide`` gives each base's parent letter; only A-T and G-C remain.
3. Watson-Crick geometry is confirmed by one of two branches:

   * branch ``H`` (every nucleic-acid base has hydrogens): both Leontis-Westhof edges
     must be WATSON_CRICK and the glycosidic bond orientation must be CIS;
   * branch ``N`` (heavy atoms only): purine N1 to pyrimidine N3 within the
     cutoff, the flanking contact (A.N6-T.O4 or G.O6-C.N4) within the
     cutoff, and CIS orientation.

4. One partner per base. Branch H lets ``base_pairs(unique=True)`` choose, because
   real H-bond geometry identifies the right partner. Without hydrogens, biotite's
   H-bond count is only a distance check and can prefer a diagonal neighbour
   (e.g. 1BF5 B1005-C2014 over B1005-C2015), so branch N asks for every candidate
   (``unique=False``) and keeps, for each base, the accepted pair with the
   shortest N1-N3 distance.

5. Pairs found so far are ``observed``. :mod:`mdpipe.infer` then adds ``inferred``
   pairs along the duplex registers they define (gap filling and outward extension,
   complementary bases only, never into a hairpin loop).
"""

import warnings

import numpy as np
import pandas as pd
import biotite.structure as struc

# The pair rule (A-T, G-C, A-U only with RNA uracil) lives in infer.is_wc_pair.
from .infer import is_wc_pair, nucleic_residues, infer_pairs

PURINES = {"A", "G"}
# Flanking Watson-Crick contact (purine atom, pyrimidine atom) per purine.
FLANK = {"A": ("N6", "O4"), "G": ("O6", "N4")}
# Shortest complementary stretch that makes another chain a likely pairing partner.
MIN_PARTNER_RUN = 4
# Columns of the accepted-pairs table; diagnostics stay in the candidates table.
PAIR_COLUMNS = ["chain1", "res_id1", "ins_code1", "res_name1", "base1",
                "chain2", "res_id2", "ins_code2", "res_name2", "base2", "pair", "evidence"]


def _atom_coord(res_atoms, name):
    """Coordinates of the named atom in a residue, or None if it is missing."""
    mask = res_atoms.atom_name == name
    return res_atoms.coord[mask][0] if mask.any() else None


def _distance(res1, name1, res2, name2):
    c1, c2 = _atom_coord(res1, name1), _atom_coord(res2, name2)
    return float(np.linalg.norm(c1 - c2)) if c1 is not None and c2 is not None else np.nan


def _base_hydrogen_count(res_atoms):
    """Number of hydrogens on the base (sugar/phosphate hydrogens are primed or HOP*)."""
    is_h = np.isin(res_atoms.element, ["H", "D"])
    on_base = np.array(["'" not in n and not n.startswith("HOP") for n in res_atoms.atom_name])
    return int((is_h & on_base).sum())


def hydrogen_status(nucleic_atoms):
    """Return ``"all"``, ``"none"`` or ``"partial"`` for hydrogens on the nucleic-acid bases."""
    counts = [_base_hydrogen_count(res) for res in struc.residue_iter(nucleic_atoms)]
    if all(c > 0 for c in counts):
        return "all"
    if not any(counts):
        return "none"
    return "partial"


def _edges_and_orientation(atoms, pairs, notes, need_edges):
    """Call biotite's edge and glycosidic-bond functions, generating bonds if needed."""
    def call():
        edges = struc.base_pairs_edge(atoms, pairs) if need_edges else None
        glyco = struc.base_pairs_glycosidic_bond(atoms, pairs)
        return edges, glyco

    try:
        return call()
    except Exception as err:  # biotite may need a BondList for the H-bond search
        if atoms.bonds is not None:
            raise
        notes.append(f"Edge/orientation call failed without bonds ({err!r}); "
                     "bonds generated from residue names and retried.")
        atoms.bonds = struc.connect_via_residue_names(atoms)
        return call()


def _complementary_partner(chain, residues, min_run=MIN_PARTNER_RUN):
    """Other chain with the longest antiparallel Watson-Crick-complementary stretch.

    Compares the sequence of ``chain`` with every other chain in ``residues`` (output of
    :func:`nucleic_residues`): residue i of one strand against residue j - k of the other
    while i + k advances, as in an antiparallel duplex.

    Returns
    -------
    (partner, length)
        The partner chain and the stretch length, or ``(None, 0)`` if no other chain has
        a complementary stretch of at least ``min_run`` nucleotides.
    """
    own = residues.get(chain, {"residues": []})["residues"]
    best, best_run = None, 0
    for other, data in residues.items():
        if other == chain:
            continue
        theirs = data["residues"]
        # Try every starting alignment and extend while consecutive bases stay complementary.
        for i in range(len(own)):
            for j in range(len(theirs)):
                run = 0
                while (i + run < len(own) and j - run >= 0
                       and is_wc_pair(own[i + run]["letter"], own[i + run]["is_rna"],
                                      theirs[j - run]["letter"], theirs[j - run]["is_rna"])):
                    run += 1
                if run > best_run:
                    best, best_run = other, run
    return (best, best_run) if best_run >= min_run else (None, 0)


def _keep_one_partner_per_base(candidates):
    """Branch N: among accepted candidates, keep each base's best pair (shortest N1-N3).

    Greedy: accepted pairs are visited from shortest to longest N1-N3; a pair is kept
    only if neither base is already used. Losers are marked in ``reject_reason``.
    Modifies ``candidates`` in place.
    """
    accepted = candidates.index[candidates["reject_reason"] == ""]
    used = {}
    for k in candidates.loc[accepted].sort_values(["n1_n3", "flank"]).index:
        row = candidates.loc[k]
        b1 = (row["chain1"], row["res_id1"], row["ins_code1"])
        b2 = (row["chain2"], row["res_id2"], row["ins_code2"])
        taken = [b for b in (b1, b2) if b in used]
        if taken:
            b = taken[0]
            candidates.loc[k, "reject_reason"] = (
                f"{b[0]}{b[1]}{b[2]} already paired with {used[b]} (shorter N1-N3)")
        else:
            used[b1] = f"{b2[0]}{b2[1]}{b2[2]}"
            used[b2] = f"{b1[0]}{b1[1]}{b1[2]}"


def find_wc_pairs(atoms, chains, cutoff=4.0):
    """Find Watson-Crick base pairs between nucleic-acid residues (DNA and/or RNA).

    Parameters
    ----------
    atoms : AtomArray
        One model, as returned by :func:`mdpipe.load_structure`.
    chains : list of str
        Nucleic-acid chain IDs (DNA, RNA or hybrid) from Part 1.
    cutoff : float
        Heavy-atom distance cutoff (Angstrom) for branch N contacts. 4.0 A keeps
        stretched pairs in medium-resolution structures (e.g. 1AOI, 1BF5 at
        2.8-2.9 A), where 3.5 A rejected about a quarter of the real pairs.

    Returns
    -------
    pairs : pandas.DataFrame
        Watson-Crick pairs, one row each: identity columns plus ``evidence``
        (``observed`` by geometry or ``inferred`` along a duplex register).
    candidates : pandas.DataFrame
        Every pair returned by ``base_pairs`` with diagnostics (C1'-C1', N1-N3 and
        flanking distances, edges, glycosidic orientation, method) and the reason it
        was rejected, if it was.
    info : dict
        ``method`` (``"H"`` or ``"N"``), ``hydrogens`` status, ``notes``, ``warnings``,
        ``n_observed`` and ``n_inferred`` pair counts, and ``registers`` (see
        :func:`mdpipe.infer.infer_pairs`).
    """
    notes, warns = [], []
    mask = struc.filter_nucleotides(atoms) & np.isin(atoms.chain_id, chains)
    nucleic = atoms[mask]
    if nucleic.array_length() == 0:
        empty = pd.DataFrame()
        return empty, empty, {"method": None, "hydrogens": None, "notes": notes, "warnings": warns,
                              "n_observed": 0, "n_inferred": 0, "registers": []}

    h_status = hydrogen_status(nucleic)
    if h_status == "partial":
        warns.append("Only some nucleic-acid bases have hydrogens; using the heavy-atom branch (N) "
                     "for the whole structure.")
    method = "H" if h_status == "all" else "N"

    raw = struc.base_pairs(nucleic, unique=(method == "H"))
    if len(raw):
        raw = np.unique(np.sort(raw, axis=1), axis=0)  # guard against (i, j) and (j, i)
    starts = struc.get_residue_starts(nucleic, add_exclusive_stop=True)
    stop_of = dict(zip(starts[:-1], starts[1:]))

    def residue(start):
        return nucleic[start:stop_of[start]]

    edges, glyco = (None, None)
    if len(raw):
        edges, glyco = _edges_and_orientation(nucleic, raw, notes, need_edges=(method == "H"))

    rows = []
    for k, (s1, s2) in enumerate(raw):
        r1, r2 = residue(s1), residue(s2)
        l1, l2 = struc.map_nucleotide(r1)[0], struc.map_nucleotide(r2)[0]
        # Report the lower (chain, residue) first.
        if (r1.chain_id[0], r1.res_id[0]) > (r2.chain_id[0], r2.res_id[0]):
            r1, r2, l1, l2 = r2, r1, l2, l1
            if edges is not None:
                edges[k] = edges[k][::-1]

        row = {
            "chain1": r1.chain_id[0], "res_id1": int(r1.res_id[0]), "ins_code1": r1.ins_code[0],
            "res_name1": r1.res_name[0], "base1": l1,
            "chain2": r2.chain_id[0], "res_id2": int(r2.res_id[0]), "ins_code2": r2.ins_code[0],
            "res_name2": r2.res_name[0], "base2": l2,
            "pair": f"{l1}-{l2}",
            "c1_c1": _distance(r1, "C1'", r2, "C1'"),
            "glycosidic": struc.GlycosidicBond(int(glyco[k])).name if glyco is not None else None,
            "edge1": struc.Edge(int(edges[k][0])).name if edges is not None else None,
            "edge2": struc.Edge(int(edges[k][1])).name if edges is not None else None,
            "n1_n3": np.nan, "flank": np.nan, "method": method, "reject_reason": "",
        }

        rna1, rna2 = bool((r1.atom_name == "O2'").any()), bool((r2.atom_name == "O2'").any())
        if not is_wc_pair(l1, rna1, l2, rna2):
            row["reject_reason"] = "not Watson-Crick (A-T/G-C/A-U)"
        else:
            pur, pyr = (r1, r2) if l1 in PURINES else (r2, r1)
            pur_letter = l1 if l1 in PURINES else l2
            row["n1_n3"] = _distance(pur, "N1", pyr, "N3")
            row["flank"] = _distance(pur, FLANK[pur_letter][0], pyr, FLANK[pur_letter][1])

            if row["glycosidic"] != "CIS":
                row["reject_reason"] = f"glycosidic bond {row['glycosidic']}"
            elif method == "H":
                if not (row["edge1"] == row["edge2"] == "WATSON_CRICK"):
                    row["reject_reason"] = f"edges {row['edge1']}/{row['edge2']}"
            else:
                if not row["n1_n3"] <= cutoff:
                    row["reject_reason"] = f"N1-N3 {row['n1_n3']:.2f} A > {cutoff}"
                elif not row["flank"] <= cutoff:
                    row["reject_reason"] = f"flanking contact {row['flank']:.2f} A > {cutoff}"
        rows.append(row)

    candidates = pd.DataFrame(rows)
    if method == "N" and len(candidates):
        _keep_one_partner_per_base(candidates)
    accepted = candidates[candidates["reject_reason"] == ""] if len(candidates) else candidates

    # Observed pairs as (chain, index) tuples, then inferred pairs along their registers.
    residues = nucleic_residues(atoms, chains)
    index_of = {(r["chain"], r["res_id"], r["ins_code"]): k
                for c in residues.values() for k, r in enumerate(c["residues"])}
    observed = [((row.chain1, index_of[(row.chain1, row.res_id1, row.ins_code1)]),
                 (row.chain2, index_of[(row.chain2, row.res_id2, row.ins_code2)]))
                for row in accepted.itertuples()]
    inferred, registers = infer_pairs(
        {c: [r["letter"] for r in v["residues"]] for c, v in residues.items()},
        {c: v["links"] for c, v in residues.items()},
        observed,
        rna={c: [r["is_rna"] for r in v["residues"]] for c, v in residues.items()})

    def pair_row(a, b, evidence):
        (c1, i), (c2, j) = sorted((a, b))
        r1, r2 = residues[c1]["residues"][i], residues[c2]["residues"][j]
        return {"chain1": c1, "res_id1": r1["res_id"], "ins_code1": r1["ins_code"],
                "res_name1": r1["res_name"], "base1": r1["letter"],
                "chain2": c2, "res_id2": r2["res_id"], "ins_code2": r2["ins_code"],
                "res_name2": r2["res_name"], "base2": r2["letter"],
                "pair": f"{r1['letter']}-{r2['letter']}", "evidence": evidence,
                "_order": (c1, i)}

    rows = ([pair_row(a, b, "observed") for a, b in observed]
            + [pair_row(a, b, "inferred") for a, b in inferred])
    pairs = (pd.DataFrame(rows).sort_values("_order")[PAIR_COLUMNS] if rows
             else pd.DataFrame(columns=PAIR_COLUMNS))

    # Chains without any pair: a note by default, a warning when another chain could pair with them.
    for chain in chains:
        if len(pairs) and ((pairs["chain1"] == chain) | (pairs["chain2"] == chain)).any():
            continue
        partner, run = _complementary_partner(chain, residues)
        if partner is None:
            notes.append(f"Nucleic-acid chain {chain} has no Watson-Crick pairs (single-stranded, "
                         "or its partner is not in the file).")
        else:
            warns.append(f"Nucleic-acid chain {chain} has no Watson-Crick pairs, although chain "
                         f"{partner} contains a complementary stretch of {run} nucleotides; "
                         "check the structure.")
    for w in warns:
        warnings.warn(w)

    return pairs.reset_index(drop=True), candidates, {
        "method": method, "hydrogens": h_status, "notes": notes, "warnings": warns,
        "n_observed": len(observed), "n_inferred": len(inferred), "registers": registers}
