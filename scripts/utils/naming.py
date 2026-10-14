"""Translate force-field residue names (AMBER, CHARMM) to PDB residue names.

biotite's residue filters and templates rely on PDB Chemical Component Dictionary
names, so files from tleap/ParmEd (DA5, DC3, RU, HID, CYX, NALA, ...) must be
translated before analysis.

* Nucleotides are identified from their atoms, independent of the name: the base
  from ``map_nucleotide``, DNA vs RNA from the O2' atom (DA5 -> DA, A3 -> A, RU -> U).
* Amino acids use a table of protonation/terminal aliases (HID -> HIS, NALA -> ALA).
"""

import numpy as np
import biotite.structure as struc

SUGAR_ATOMS = {"C1'", "C4'", "O4'"}
PDB_NUCLEOTIDES = {"DA", "DC", "DG", "DT", "DU", "A", "C", "G", "U"}

STANDARD_AMINO_ACIDS = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
AMINO_ACID_ALIASES = {
    # AMBER protonation / disulfide states
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HIN": "HIS",
    "CYX": "CYS", "CYM": "CYS", "ASH": "ASP", "GLH": "GLU", "LYN": "LYS",
    # CHARMM histidines
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
}


def amino_acid_name(name):
    """PDB name for an amino-acid residue name, or the name unchanged if not recognised.

    Handles protonation aliases and AMBER terminal prefixes (NALA, CHID -> ALA, HIS).
    """
    if name in STANDARD_AMINO_ACIDS:
        return name
    if name in AMINO_ACID_ALIASES:
        return AMINO_ACID_ALIASES[name]
    if len(name) == 4 and name[0] in "NC":
        core = AMINO_ACID_ALIASES.get(name[1:], name[1:])
        if core in STANDARD_AMINO_ACIDS:
            return core
    return name


def is_forcefield_nucleotide_name(name):
    """True for standard or force-field nucleotide names: A, DA, RA, DA5, A3, RU, DAN, ...

    Modified nucleotides with their own dictionary names (PSU, 5MC, H2U, ...) are not
    matched, so they keep their identity and their own templates.
    """
    core = name[1:] if name[:1] in "DR" and len(name) > 1 else name
    core = core[:-1] if core[-1:] in "53N" and len(core) > 1 else core
    return core in {"A", "C", "G", "T", "U"}


def pdb_residue_name(res):
    """PDB name for one residue (AtomArray of a single residue).

    Standard/force-field nucleotide names -> DA/DC/DG/DT/DU or A/C/G/U, decided from the
    atoms; amino acids via :func:`amino_acid_name`; anything else (modified residues,
    ligands) keeps its name.
    """
    names = set(res.atom_name)
    if SUGAR_ATOMS <= names and is_forcefield_nucleotide_name(res.res_name[0]):
        letter = struc.map_nucleotide(res)[0]
        if letter is not None:
            return letter if "O2'" in names else "D" + letter
        return res.res_name[0]
    return amino_acid_name(res.res_name[0])


def translate_residue_names(atoms):
    """Copy of ``atoms`` with every residue renamed by :func:`pdb_residue_name`."""
    atoms = atoms.copy()
    starts = struc.get_residue_starts(atoms, add_exclusive_stop=True)
    names = atoms.res_name.astype(object)
    for a, b in zip(starts[:-1], starts[1:]):
        names[a:b] = pdb_residue_name(atoms[a:b])
    atoms.res_name = np.array(names, dtype=atoms.res_name.dtype)
    return atoms
