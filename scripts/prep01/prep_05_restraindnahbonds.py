#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Add distance restraints on Watson-Crick H-bonds of DNA/RNA base pairs.

Base pairs are detected from the coordinates (``.gro``) rather than assumed from the
sequences, so overhangs, hairpins, RNA and duplexes split over any number of
molecules are handled without manual residue ranges:

1. ``[ molecules ]`` + the included ``.itp`` files give every molecule instance and its
   global atom range (GROMACS numbers atoms in ``[ molecules ]`` order).
2. The ``.gro`` is checked atom by atom against the topology up to the last
   nucleic-acid molecule (the solvent/ion tail may differ, e.g. a pre-genion ``.gro``).
3. The nucleotides of those molecules are grouped into strands (runs joined by O3'-P
   links); amino acids, ligands and isolated nucleotides are left out. Residue names
   are translated to PDB names, and ``scripts.utils.nucleic`` finds the pairs
   (observed by geometry + inferred along duplex registers). A system without any
   Watson-Crick pair (e.g. a single strand) needs no restraints and is not an error.
4. Restraints between atoms of different molecules go to
   ``[ intermolecular_interactions ]`` in ``topol`` (global numbering); restraints within
   one molecule (hairpins, merged strands) go to that molecule's ``.itp`` (local
   numbering, applied to every copy of the molecule).
"""

import argparse
import io
import itertools
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import biotite.structure as struc
import biotite.structure.io.gro as gro_io

from scripts.utils.helpers import setup_logger
from scripts.utils.nucleic import classify_structure, nucleic_chain_ids, find_wc_pairs, pairs_diagram
from scripts.utils.nucleic.infer import LINK_MAX
from scripts.utils.naming import SUGAR_ATOMS, pdb_residue_name

# H-bonded atoms per Watson-Crick pair, as (purine atom, pyrimidine atom).
# A pairs with T or (RNA) U; both pyrimidines use O4 and N3.
RESTRAINT_ATOMS = {"A": [("N6", "O4"), ("N1", "N3")],
                   "G": [("N2", "O2"), ("N1", "N3"), ("O6", "N4")]}
CHECKPOINT = "check_05_restraindnahbonds_*.done"
# Shortest run of linked nucleotides treated as a strand that can pair.
MIN_STRAND_LENGTH = 2
# Columns of the restraint review table (restraints_hbonds.csv).
RESTRAINT_COLUMNS = ["label", "pair", "evidence", "chain1", "res_id1", "res_name1", "chain2",
                     "res_id2", "res_name2", "atom_pur", "atom_pyr", "global_ai", "global_aj",
                     "section"]


def main(logfile: Path, gmxdir: Path, topol: Path = Path("./topol.top"), gro: Path = None,
         type: int = 1, type2: int = 2, low: float = 0.0, up1: float = 0.3, up2: float = 0.4,
         fac: float = 1.0):
    """
    Detect Watson-Crick base pairs (DNA and RNA) and add H-bond distance restraints.

    Arguments:
        logfile : Path to the log file.
        gmxdir  : Directory containing GROMACS files (with toppar/).
        topol   : Topology file to append restraints to (default: ./topol.top).
        gro     : Coordinates in topology atom order (default: the single *.gro in gmxdir).
        type    : Restraint type (default: 1)
        type2   : Secondary restraint type (default: 2)
        low     : Lower bound distance (nm) (default: 0.0)
        up1     : Upper bound distance 1 (nm) (default: 0.3)
        up2     : Upper bound distance 2 (nm) (default: 0.4)
        fac     : Force constant scaling factor (default: 1.0)
    """
    logger = setup_logger("prep", logfile)
    gmxdir, topol = Path(gmxdir), Path(topol)
    toppar = gmxdir / "toppar"

    ##### CHECK INPUT #####
    if not gmxdir.is_dir():
        logger.error(f"Directory {gmxdir} not found")
        sys.exit(1)
    if not toppar.is_dir():
        logger.error(f"Topology directory {toppar} not found")
        sys.exit(1)
    if not topol.is_file():
        logger.error(f"Topology file {topol} not found.")
        sys.exit(1)
    if low > up1 or up1 > up2:
        logger.error("Invalid distance bounds: ensure that low <= up1 <= up2.")
        sys.exit(1)
    if fac <= 0:
        logger.error("Force constant scaling factor (fac) should be positive.")
        sys.exit(1)
    if gro is None:
        candidates = sorted(gmxdir.glob("*.gro"))
        if len(candidates) != 1:
            logger.error(f"Expected exactly one .gro in {gmxdir} (found {len(candidates)}); "
                         "pass gro explicitly.")
            sys.exit(1)
        gro = candidates[0]
    gro = Path(gro)
    if not gro.is_file():
        logger.error(f"Coordinate file {gro} not found.")
        sys.exit(1)

    ##### READ TOPOLOGY #####
    includes, molecules = _read_topol(topol, logger)
    itp_paths = [_resolve_include(inc, topol, toppar) for inc in includes]
    moltypes = {}
    for path in itp_paths:
        if path is None:
            continue
        for name, atoms in _read_moleculetypes(path).items():
            moltypes[name] = {"itp": path, "atoms": atoms}

    # Checkpoint: restraints may already be in topol or in any molecule itp.
    if list(gmxdir.glob(CHECKPOINT)):
        logger.debug("Checkpoint file found.")
        written = [p for p in [topol] + [m["itp"] for m in moltypes.values()]
                   if "[ distance_restraints ]" in p.read_text()]
        if written:
            logger.info(f"[ distance_restraints ] found in: {', '.join(str(p) for p in written)}.")
            logger.warning("Skipping H-bond restraints addition.")
            return
        logger.warning("Checkpoint file found but no [ distance_restraints ] section; the previous "
                       "run might have been incomplete. Proceeding.")

    # Global atom ranges per molecule instance, in [ molecules ] order.
    instances, start = [], 0
    for name, count in molecules:
        if name not in moltypes:
            logger.error(f"Molecule {name} from [ molecules ] not defined in any included itp.")
            sys.exit(1)
        natoms = len(moltypes[name]["atoms"])
        is_nucleic = _is_nucleic(moltypes[name]["atoms"])
        for copy in range(count):
            instances.append({"moltype": name, "copy": copy, "start": start, "natoms": natoms,
                              "nucleic": is_nucleic})
            start += natoms
    nucleic = [m for m in instances if m["nucleic"]]
    if not nucleic:
        logger.warning("Nucleic-acid molecules are expected but could not be found in the topology. Hydrogen bond restrains could not be added. This may cause errors downstream.")
        return
    logger.debug("Nucleic-acid molecules: " + ", ".join(
        f"{m['moltype']}#{m['copy'] + 1} (atoms {m['start'] + 1}-{m['start'] + m['natoms']})"
        for m in nucleic))

    ##### READ AND CHECK COORDINATES #####
    # Only atoms up to the end of the last nucleic-acid molecule are needed.
    n_check = nucleic[-1]["start"] + nucleic[-1]["natoms"]
    atoms, gro_total = _read_gro_prefix(gro, n_check)
    if gro_total < n_check:
        logger.error(f"{gro} has {gro_total} atoms; the topology needs at least "
                     f"{n_check} to cover the nucleic acids.")
        sys.exit(1)
    # Expected names only for the molecules up to the last nucleic acid (not the solvent).
    topo_names = [(a["residue"], a["atom"]) for m in instances if m["start"] < n_check
                  for a in moltypes[m["moltype"]]["atoms"]]
    for k, (res, atom) in enumerate(topo_names):
        if (atoms.res_name[k], atoms.atom_name[k]) != (res[:5], atom[:5]):
            logger.error(f"Atom {k + 1}: {gro.name} has {atoms.res_name[k]} {atoms.atom_name[k]}, "
                         f"topology has {res} {atom}. Coordinates and topology do not match.")
            sys.exit(1)
    # The total atom count comes from the .gro header, so the solvent itself is never read.
    total = start
    if gro_total != total:
        logger.warning(f"{gro.name} has {gro_total} atoms, the topology {total}.")

    ##### NUCLEIC-ACID STRANDS #####
    # Keep only nucleotides, grouped into strands with one chain label each.
    try:
        sub, global_index, strand_info, ignored = _nucleic_strands(atoms, nucleic)
    except RuntimeError as err:
        logger.error(str(err))
        sys.exit(1)
    if ignored:
        logger.debug("Ignored single nucleotides (not part of a strand): " + ", ".join(ignored))
    if sub is None:
        logger.error("No nucleic-acid strands (at least two linked nucleotides) found, but were expected. Hydrogen bond restrains could not be added. This may cause errors downstream.")
        return
    for chain, (mol, n_res) in strand_info.items():
        logger.debug(f"Strand {chain}: {mol['moltype']}#{mol['copy'] + 1}, {n_res} nucleotides")

    ##### FIND BASE PAIRS #####
    chains_df, summary = classify_structure(sub)
    nuc_chains = nucleic_chain_ids(chains_df)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pairs, candidates, info = find_wc_pairs(sub, nuc_chains)
    # Per-strand notes (e.g. unpaired strands) go to INFO, likely problems to WARNING.
    for note in info["notes"]:
        logger.debug(note)
    for w in info["warnings"]:
        logger.warning(w)
    logger.debug(f"Nucleic acids: {summary}; branch {info['method']} (hydrogens: {info['hydrogens']}); "
                f"{len(pairs)} Watson-Crick pairs ({info['n_observed']} observed, "
                f"{info['n_inferred']} inferred).")
    if not len(pairs):
        logger.debug("No Watson-Crick base pairs found. No restraints needed.")

    ##### RESTRAINT LINES #####
    atom_index = {(sub.chain_id[k], int(sub.res_id[k]), sub.atom_name[k]): global_index[k]
                  for k in range(sub.array_length())}
    inter, intra = [], {}     # intra: moltype -> lines
    rows = []
    for label, p in enumerate(pairs.itertuples(), start=1):
        # Describes each base of the pair as (strand label, residue number, base letter), e.g. ("A", 1155, "G")
        side1 = (p.chain1, int(p.res_id1), p.base1)
        side2 = (p.chain2, int(p.res_id2), p.base2)
        # The purine is the one with a defined restraint atom; the pyrimidine is the other.
        # Ordering matters because the atom table is written as (purine atom, pyrimidine atom)
        pur, pyr = (side1, side2) if p.base1 in RESTRAINT_ATOMS else (side2, side1)
        # Check which molecule copy each strand came from
        mol1, mol2 = strand_info[p.chain1][0], strand_info[p.chain2][0]
        # Loop over the H-bond atom pairs for that purine
        for pur_atom, pyr_atom in RESTRAINT_ATOMS[pur[2]]:
            # Convert each atom, identified by (strand, residue, atom name), into its global atom index
            ai = atom_index.get((pur[0], pur[1], pur_atom))
            aj = atom_index.get((pyr[0], pyr[1], pyr_atom))
            if ai is None or aj is None:
                logger.warning(f"Missing atom for pair {p.chain1}{p.res_id1}-{p.chain2}{p.res_id2} "
                             f"({pur_atom}/{pyr_atom}).")
                continue
            # If the two strands are from the same molecule, use local numbering in its itp; otherwise, use global numbering in topol.
            if mol1 is mol2:     
                off = mol1["start"]
                intra.setdefault(mol1["moltype"], []).append(
                    _restraint_line(ai - off + 1, aj - off + 1, type, label, type2, low, up1, up2, fac))
            else:
                inter.append(_restraint_line(ai + 1, aj + 1, type, label, type2, low, up1, up2, fac))
            rows.append({"label": label, "pair": p.pair, "evidence": p.evidence,
                         "chain1": p.chain1, "res_id1": p.res_id1, "res_name1": p.res_name1,
                         "chain2": p.chain2, "res_id2": p.res_id2, "res_name2": p.res_name2,
                         "atom_pur": pur_atom, "atom_pyr": pyr_atom,
                         "global_ai": ai + 1, "global_aj": aj + 1,
                         "section": "intramolecular" if mol1 is mol2 else "intermolecular"})

    ##### WRITE (only after everything above succeeded) #####
    header = "; ai    aj    type    index    type'    low    up1    up2    fac\n"
    if inter:
        with open(topol, "a") as f:
            f.write("\n[ intermolecular_interactions ]\n[ distance_restraints ]\n" + header)
            f.writelines(inter)
        logger.info(f"{len(inter)} intermolecular H-bond restraints added to {topol}.")
    for name, lines in intra.items():
        itp = moltypes[name]["itp"]
        count = dict(molecules)[name]
        if count > 1:
            logger.debug(f"{name} appears {count} times; its intramolecular restraints apply to every copy.")
        with open(itp, "a") as f:
            f.write("\n[ distance_restraints ]\n" + header)
            f.writelines(lines)
        logger.debug(f"{len(lines)} intramolecular H-bond restraints added to {itp} ({name}).")

    # Review files: pair diagram and restraint table.
    pd.DataFrame(rows, columns=RESTRAINT_COLUMNS).to_csv(gmxdir / "restraints_hbonds.csv",
                                                          index=False)
    legend = "\n".join(f"- strand {c}: {m['moltype']} copy {m['copy'] + 1}, {n} nucleotides"
                       for c, (m, n) in strand_info.items())
    (gmxdir / "restraints_basepairs.md").write_text(
        pairs_diagram(sub, nuc_chains, pairs, info, title=gro.stem)
        + "\n## Strands\n\n" + legend + "\n", encoding="utf-8")
    logger.debug(f"Files written: {gmxdir / 'restraints_basepairs.md'}, "
                f"{gmxdir / 'restraints_hbonds.csv'}")

    checkpoint = gmxdir / f"check_05_restraindnahbonds_{int(time.time())}.done"
    checkpoint.write_text(f"H-bond restraint step finished: {len(rows)} restraint(s) for "
                          f"{len(pairs)} Watson-Crick pair(s).\n")


def _read_gro_prefix(gro, n_atoms):
    """First ``n_atoms`` atoms of a .gro file, and the total atom count from its header.

    Only the title, the atom-count line and the first ``n_atoms`` atom lines are read;
    the rest of the file (typically solvent and ions) is skipped.

    Returns
    -------
    atoms : AtomArray
        The first ``min(n_atoms, total)`` atoms.
    total : int
        Number of atoms stated on line 2 of the file.
    """
    with open(gro) as f:
        title = f.readline()
        total = int(f.readline().split()[0])
        # One atom per line, in file order: stop after the atoms that are needed.
        lines = list(itertools.islice(f, min(n_atoms, total)))
    # Rebuild a small, valid .gro (the box line is a placeholder; it is not used here).
    text = title + f"{len(lines)}\n" + "".join(lines) + "   0.00000   0.00000   0.00000\n"
    atoms = gro_io.GROFile.read(io.StringIO(text)).get_structure(model=1)
    return atoms, total


def _restraint_line(ai, aj, type, label, type2, low, up1, up2, fac):
    return f"{ai}    {aj}    {type}    {label}    {type2}    {low}    {up1}    {up2}    {fac}\n"


def _read_topol(topol, logger):
    """Included files (in order) and the [ molecules ] list of (name, count)."""
    includes, molecules, section = [], [], None
    for raw in open(topol):
        line = raw.split(";")[0].strip()
        if not line:
            continue
        if line.startswith("#include"):
            includes.append(line.split(None, 1)[1].strip().strip('"'))
        elif line.startswith("["):
            section = line.strip("[] ").lower()
        elif section == "molecules":
            cols = line.split()
            if len(cols) == 2 and cols[1].isdigit():
                molecules.append((cols[0], int(cols[1])))
            else:
                logger.warning(f"[ molecules ] line '{line}' in {topol} not understood; skipped.")
    return includes, molecules


def _resolve_include(include, topol, toppar):
    """Path of an included itp (relative to topol, else in toppar); None if not found."""
    for path in (topol.parent / include, toppar / Path(include).name):
        if path.is_file():
            return path
    return None


def _read_moleculetypes(path):
    """{moleculetype name: [atom dicts]} for every moleculetype defined in an itp."""
    out, name, section, expect_name = {}, None, None, False
    for raw in open(path):
        line = raw.split(";")[0].strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("["):
            section = line.strip("[] ").lower()
            expect_name = section == "moleculetype"
            continue
        if expect_name:
            name = line.split()[0]
            out[name] = []
            expect_name = False
        elif section == "atoms" and name is not None:
            cols = line.split()
            out[name].append({"nr": int(cols[0]), "resnr": int(cols[2]),
                              "residue": cols[3], "atom": cols[4]})
    return out


def _is_nucleic(atom_list):
    """True if any residue of the moleculetype has a nucleotide sugar."""
    by_res = {}
    for a in atom_list:
        by_res.setdefault(a["resnr"], set()).add(a["atom"])
    return any(SUGAR_ATOMS <= names for names in by_res.values())


def _nucleic_strands(atoms, nucleic, min_length=MIN_STRAND_LENGTH):
    """Nucleotide strands of the nucleic-acid molecules, one chain label per strand.

    Only nucleotides are kept (sugar atoms and a recognisable base); amino acids, ligands
    and other residues inside these molecules are left out. A strand is a run of
    nucleotides joined by O3'-P links; strands shorter than ``min_length`` (isolated
    nucleotides such as nucleotide cofactors) are left out as well.

    Returns
    -------
    sub : AtomArray
        The kept nucleotides, with PDB residue names and one chain label per strand.
    global_index : ndarray
        Global (0-based) atom index of each atom of ``sub``.
    strand_info : dict
        ``{chain: (molecule instance, number of nucleotides)}``.
    ignored : list of str
        Nucleotides left out because their strand is too short.
    """
    pieces, index, strand_info, ignored = [], [], {}, []
    labels = _chain_ids()
    for mol in nucleic:
        part = atoms[mol["start"]:mol["start"] + mol["natoms"]].copy()
        part.ins_code[:] = ""
        starts = struc.get_residue_starts(part, add_exclusive_stop=True)

        # Group consecutive linked nucleotides into strands; any other residue ends a strand.
        strands, current, prev = [], [], None
        for a, b in zip(starts[:-1], starts[1:]):
            res = part[a:b]
            if not _is_nucleotide(res):
                if current:
                    strands.append(current)
                current, prev = [], None
                continue
            if prev is not None and not _linked(prev, res):
                strands.append(current)
                current = []
            current.append((a, b))
            prev = res
        if current:
            strands.append(current)

        for strand in strands:
            # Too short to be part of a duplex: record it and skip it.
            if len(strand) < min_length:
                ignored += [f"{part.res_name[a]}{part.res_id[a]} ({mol['moltype']}#{mol['copy'] + 1})"
                            for a, _ in strand]
                continue
            chain = next(labels, None)
            if chain is None:
                raise RuntimeError("More nucleic-acid strands than available chain labels.")

            # Copy the strand's atoms, label them and translate residue names to PDB names.
            sel = np.concatenate([np.arange(a, b) for a, b in strand])
            piece = part[sel]
            piece.chain_id[:] = chain
            res_starts = struc.get_residue_starts(piece, add_exclusive_stop=True)
            for a, b in zip(res_starts[:-1], res_starts[1:]):
                piece.res_name[a:b] = pdb_residue_name(piece[a:b])
            pieces.append(piece)
            index.append(mol["start"] + sel)
            strand_info[chain] = (mol, len(strand))

    if not pieces:
        return None, np.array([], dtype=int), strand_info, ignored
    sub = struc.concatenate(pieces) if len(pieces) > 1 else pieces[0]
    return sub, np.concatenate(index), strand_info, ignored


def _is_nucleotide(res):
    """True if the residue has nucleotide sugar atoms and a base biotite can identify."""
    return SUGAR_ATOMS <= set(res.atom_name) and struc.map_nucleotide(res)[0] is not None


def _linked(res_a, res_b):
    """True if res_a's O3' is covalently linked to res_b's P (distance <= LINK_MAX)."""
    o3, p = res_a.coord[res_a.atom_name == "O3'"], res_b.coord[res_b.atom_name == "P"]
    return len(o3) > 0 and len(p) > 0 and np.linalg.norm(o3[0] - p[0]) <= LINK_MAX


def _chain_ids():
    alphabet = [chr(c) for c in range(ord("A"), ord("Z") + 1)]
    # One-, two- and three-letter labels, in that order.
    for length in (1, 2, 3):
        for letters in itertools.product(alphabet, repeat=length):
            yield "".join(letters)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--logfile", type=Path, required=True)
    parser.add_argument("--gmxdir", type=Path, required=True)
    parser.add_argument("--topol", type=Path, default=Path("./topol.top"))
    parser.add_argument("--gro", type=Path, default=None)
    parser.add_argument("--type", type=int, default=1)
    parser.add_argument("--type2", type=int, default=2)
    parser.add_argument("--low", type=float, default=0.0)
    parser.add_argument("--up1", type=float, default=0.3)
    parser.add_argument("--up2", type=float, default=0.4)
    parser.add_argument("--fac", type=float, default=1.0)
    main(**vars(parser.parse_args()))
