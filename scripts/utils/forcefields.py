"""Choice of AMBER force fields (LEaP leaprc names) from the molecule types present.

Each name is sourced in LEaP as ``source leaprc.<name>``. The water model is always
included and loaded last, after the biomolecular force fields.
"""

# Force field for each molecule type, in loading order.
FORCEFIELDS = {
    "protein": "protein.ff19SB",
    "phos": "phosaa19SB",          # phosphorylated amino acids; builds on ff19SB
    "dna": "DNA.bsc1",
    "rna": "RNA.OL3",
}
# Water model; the solvation step (OPCBOX, OPC ion parameters) depends on it.
WATER = "water.opc"

# Name prefixes that identify each molecule type in a user-provided list.
PREFIXES = {
    "protein": "protein.",
    "phos": "phosaa",
    "dna": "DNA.",
    "rna": "RNA.",
    "water": "water.",
}


def select_forcefields(protein: bool = False, phos: bool = False, dna: bool = False,
                       rna: bool = False, ffs: list | None = None, logger=None) -> list:
    """Force fields to source in LEaP for the given molecule types.

    Without ``ffs``, one force field per present type is chosen from FORCEFIELDS,
    followed by the water model. With ``ffs``, that list is used as given; the water
    model is appended if the list has no water entry, and a warning is logged for
    each present molecule type without a matching entry.
    """
    present = {"protein": protein, "phos": phos, "dna": dna, "rna": rna}

    # Automatic choice: one force field per present type, water last.
    if not ffs:
        return [FORCEFIELDS[kind] for kind, flag in present.items() if flag] + [WATER]

    chosen = list(ffs)
    # Detected molecule types that the given list does not cover.
    for kind, flag in present.items():
        if flag and not any(ff.startswith(PREFIXES[kind]) for ff in chosen):
            if logger is not None:
                logger.warning(
                    f"The structure contains {kind} residues, but no '{PREFIXES[kind]}*' "
                    f"force field is in ffs {chosen}. LEaP is likely to fail or leave "
                    f"atoms untyped (default for {kind}: {FORCEFIELDS[kind]}).")
    # The water model is always loaded, after the given force fields.
    if not any(ff.startswith(PREFIXES["water"]) for ff in chosen):
        chosen.append(WATER)
    return chosen
