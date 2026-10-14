"""Composition of a structure: which polymer chains it contains and of what kind.

The record lists every protein/nucleic-acid chain of the input PDB with its kind, its
position in order, residue counts per type, sequence and phosphorylated residues. The
system-level flags (``has_protein``, ``has_dna``, ``has_rna``, ``has_phospho``) are
independent and computed from residue counts, so any combination is represented.

The record is saved as JSON for review and reloaded to detect changes of the input
between runs.
"""

import hashlib
import json
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import NamedTuple

import numpy as np
import biotite.structure as struc
from biotite.structure.info import one_letter_code

from .naming import translate_residue_names
from .nucleic.classify import classify_structure, residue_type
from .nucleic.io import load_structure

SCHEMA_VERSION = 1
# File name of the record, written in the mutation directory next to the input PDB.
RECORD_NAME = "system_record.json"


@dataclass(frozen=True)
class ChainInfo:
    """One polymer chain of the input structure."""
    id: str                     # chain ID in the input PDB
    order: int                  # 1-based position among the polymer chains
    kind: str                   # protein, DNA, RNA or DNA-RNA hybrid (majority of residues)
    n_residues: int
    n_protein: int              # residues of each type in this chain
    n_dna: int
    n_rna: int
    sequence: str               # one letter per residue, X if unknown
    phospho_residues: tuple = ()  # phosphorylated amino acids, e.g. ("PTR701",)
    first_res_id: int | None = None  # residue number of the chain's first residue


@dataclass(frozen=True)
class Composition:
    """Composition of a structure file."""
    source: str                 # path of the analysed PDB
    source_sha256: str          # fingerprint of that file
    detected_at: str            # ISO timestamp of the detection
    chains: tuple = field(default_factory=tuple)
    schema_version: int = SCHEMA_VERSION
    # Input chain -> simulation chain mapping, added by the analysis stage (or None).
    chain_mapping: dict | None = None
    # Analysis settings (parameter file with resolved flags), added by the analysis stage.
    params: dict | None = None

    # System-level flags, computed from residue counts so mixed chains are not missed.
    @property
    def has_protein(self) -> bool:
        return any(c.n_protein for c in self.chains)

    @property
    def has_dna(self) -> bool:
        return any(c.n_dna for c in self.chains)

    @property
    def has_rna(self) -> bool:
        return any(c.n_rna for c in self.chains)

    @property
    def has_phospho(self) -> bool:
        return any(c.phospho_residues for c in self.chains)

    @property
    def nucleic_chains(self) -> list:
        """IDs of chains that contain DNA or RNA residues."""
        return [c.id for c in self.chains if c.n_dna or c.n_rna]

    def summary(self) -> str:
        """One line, e.g. 'protein x2 (phosphorylated: A:PTR701), DNA x2'."""
        # Number of chains of each kind, in order of first appearance.
        counts = {}
        for c in self.chains:
            counts[c.kind] = counts.get(c.kind, 0) + 1
        phospho = [f"{c.id}:{r}" for c in self.chains for r in c.phospho_residues]
        parts = []
        for kind, n in counts.items():
            text = f"{kind} x{n}"
            # Phosphorylated residues are listed next to the protein chains.
            if kind == "protein" and phospho:
                text += f" (phosphorylated: {', '.join(phospho)})"
            parts.append(text)
        return ", ".join(parts) if parts else "no protein or nucleic-acid chains"


class _Residue(NamedTuple):
    """One polymer residue of a structure, in file order."""
    chain: str
    res_id: int
    label: str          # residue name, number and insertion code, e.g. "PTR701"
    letter: str         # one-letter code, X if unknown
    kind: str           # protein, DNA or RNA
    phospho: bool       # amino acid carrying a phosphorus atom
    first_atom: int     # 0-based index of the residue's first atom in the structure
    last_atom: int      # 0-based index of its last atom


def _polymer_atoms(pdbpath):
    """First model without solvent and ions, with standard residue names.

    Returns the atoms and, for each of them, its 0-based index in the structure as read.
    """
    atoms, _ = load_structure(pdbpath, model=1)
    # Solvent and monoatomic ions are not polymer chains; dropping them keeps large
    # solvated structures fast to analyse.
    index = np.nonzero(~(struc.filter_solvent(atoms) | struc.filter_monoatomic_ions(atoms)))[0]
    # Standard residue names, so force-field names (DA5, HIE, ...) are recognised.
    return translate_residue_names(atoms[index]), index


def _polymer_residues(atoms, index) -> list:
    """Protein and nucleic-acid residues of ``atoms`` in file order (ligands skipped)."""
    residues = []
    starts = struc.get_residue_starts(atoms, add_exclusive_stop=True)
    for a, b in zip(starts[:-1], starts[1:]):
        res = atoms[a:b]
        kind = residue_type(res.res_name[0], res)
        if kind not in ("protein", "DNA", "RNA"):
            continue                                      # ligands and other molecules
        residues.append(_Residue(
            chain=str(res.chain_id[0]), res_id=int(res.res_id[0]),
            label=f"{res.res_name[0]}{res.res_id[0]}{res.ins_code[0]}",
            letter=one_letter_code(res.res_name[0]) or "X", kind=kind,
            # An amino acid carrying a phosphorus atom is phosphorylated.
            phospho=bool(kind == "protein" and (res.element == "P").any()),
            first_atom=int(index[a]), last_atom=int(index[b - 1])))
    return residues


def detect(pdbpath) -> Composition:
    """Detect the composition of a PDB file (first model)."""
    pdbpath = Path(pdbpath)
    atoms, index = _polymer_atoms(pdbpath)
    labels = dict(classify_structure(atoms)[0][["chain", "label"]].values)

    # Polymer residues grouped by chain ID, chains in order of first appearance.
    by_chain = {}
    for res in _polymer_residues(atoms, index):
        by_chain.setdefault(res.chain, []).append(res)

    chains = []
    for chain_id, residues in by_chain.items():
        # Residue counts per type, so a chain mixing types is fully described.
        counts = {kind: sum(r.kind == kind for r in residues)
                  for kind in ("protein", "DNA", "RNA")}
        chains.append(ChainInfo(
            id=chain_id, order=len(chains) + 1, kind=labels.get(chain_id, "other"),
            n_residues=len(residues), n_protein=counts["protein"],
            n_dna=counts["DNA"], n_rna=counts["RNA"],
            sequence="".join(r.letter for r in residues),
            phospho_residues=tuple(r.label for r in residues if r.phospho),
            first_res_id=residues[0].res_id))

    return Composition(source=str(pdbpath), source_sha256=file_sha256(pdbpath),
                       detected_at=datetime.now().isoformat(timespec="seconds"),
                       chains=tuple(chains))


def save(comp: Composition, path) -> None:
    """Write the record as indented JSON, with the derived flags for readability."""
    data = asdict(comp)
    data["flags"] = {"has_protein": comp.has_protein, "has_dna": comp.has_dna,
                     "has_rna": comp.has_rna, "has_phospho": comp.has_phospho}
    # Values JSON has no type for (dates, paths) are written as text.
    Path(path).write_text(json.dumps(data, indent=2, default=str) + "\n", encoding="utf-8")


def load(path) -> Composition:
    """Read a record written by :func:`save` (derived flags are recomputed, not read)."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    data.pop("flags", None)
    chains = tuple(ChainInfo(**{**c, "phospho_residues": tuple(c["phospho_residues"])})
                   for c in data.pop("chains"))
    return Composition(chains=chains, **data)


def map_to_structure(old: Composition, pdbpath) -> dict:
    """Match the polymer chains of a structure of the same system to a record, by order.

    The polymer residues of ``pdbpath`` (first model) are split into chains in file
    order: one chain per chain ID, and a run of residues without chain ID is cut into
    pieces with the lengths of the record's chains at the matching positions. The k-th
    chain is paired with the k-th chain of ``old`` and verified by length and sequence.

    Returns ``{"0": entry, "1": entry, ...}`` (0-based position), where each entry has
    ``id_old``, ``id_new``, ``offset`` (original residue number = new + offset),
    ``verified``, ``resids_new`` and ``atoms_new`` (first and last residue number and
    0-based atom index in the structure) and ``note`` (why a chain is not verified).
    Unverified chains have ``id_old`` and ``offset`` set to None.
    """
    pdbpath = Path(pdbpath)
    atoms, index = _polymer_atoms(pdbpath)
    residues = _polymer_residues(atoms, index)

    # Consecutive residues with the same chain ID form one group.
    groups = []
    for res in residues:
        if groups and groups[-1][0] == res.chain:
            groups[-1][1].append(res)
        else:
            groups.append((res.chain, [res]))

    # Each named group is one chain; a group without chain ID is cut by the lengths
    # of the record's chains at the matching positions.
    segments = []
    for chain_id, group in groups:
        if chain_id.strip():
            segments.append((chain_id, group))
            continue
        while group:
            k = len(segments)
            size = old.chains[k].n_residues if k < len(old.chains) else len(group)
            segments.append((chain_id, group[:size]))
            group = group[size:]

    mapping = {}
    for k in range(max(len(old.chains), len(segments))):
        ref = old.chains[k] if k < len(old.chains) else None
        seg = segments[k][1] if k < len(segments) else None
        entry = {"id_old": None, "id_new": None, "offset": None, "verified": False,
                 "resids_new": None, "atoms_new": None, "note": ""}
        # Location of the chain in the structure, for selection without a chain ID.
        if seg is not None:
            entry.update(id_new=segments[k][0],
                         resids_new=[seg[0].res_id, seg[-1].res_id],
                         atoms_new=[seg[0].first_atom, seg[-1].last_atom])
        # The pair is verified only if both chains exist with the same residues.
        if ref is None:
            entry["note"] = f"no chain at position {k} in {old.source}"
        elif seg is None:
            entry["note"] = f"input chain {ref.id} has no counterpart in {pdbpath.name}"
        elif len(seg) != ref.n_residues:
            entry["note"] = (f"input chain {ref.id} has {ref.n_residues} residues, "
                             f"this chain {len(seg)}")
        elif "".join(r.letter for r in seg) != ref.sequence:
            entry["note"] = f"sequence differs from input chain {ref.id}"
        elif ref.first_res_id is None:
            entry["note"] = f"no residue numbers recorded for input chain {ref.id}"
        else:
            entry.update(id_old=ref.id, offset=ref.first_res_id - seg[0].res_id,
                         verified=True)
        mapping[str(k)] = entry
    return mapping


def differences(old: Composition, new: Composition) -> list:
    """Readable list of what differs between two records (empty if equivalent)."""
    out = []
    if old.source_sha256 != new.source_sha256:
        out.append(f"input file content changed ({old.source} -> {new.source})")
    old_chains = {c.id: c for c in old.chains}
    new_chains = {c.id: c for c in new.chains}
    # Chains present in only one of the records.
    for cid in old_chains.keys() - new_chains.keys():
        out.append(f"chain {cid} no longer present")
    for cid in new_chains.keys() - old_chains.keys():
        out.append(f"chain {cid} is new")
    # Chains present in both: compare their recorded properties.
    for cid in sorted(old_chains.keys() & new_chains.keys()):
        a, b = asdict(old_chains[cid]), asdict(new_chains[cid])
        for key in a:
            if a[key] != b[key]:
                shown = "differs" if key == "sequence" else f"{a[key]} -> {b[key]}"
                out.append(f"chain {cid}: {key} {shown}")
    return out


def file_sha256(path) -> str:
    """SHA-256 of a file, read in blocks."""
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()
