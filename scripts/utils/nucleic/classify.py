"""Part 1: classify chains and the whole structure as DNA, RNA, both or neither."""

import pandas as pd
import biotite.structure as struc
from biotite.structure.info import link_type


def residue_type(res_name, res_atoms):
    """Return ``"DNA"``, ``"RNA"``, ``"protein"`` or ``None`` for one residue.

    The PDB Chemical Component Dictionary link type is used first, which covers
    modified residues. Unknown nucleotides fall back to sugar chemistry:
    an O2' atom means RNA, a C2' without O2' means DNA.
    """
    ltype = link_type(res_name)
    if ltype is not None:
        ltype = ltype.upper()
        if "DNA" in ltype:
            return "DNA"
        if "RNA" in ltype:
            return "RNA"
        if "PEPTIDE" in ltype:
            return "protein"

    names = set(res_atoms.atom_name)
    if "C1'" in names and "C4'" in names:
        if "O2'" in names:
            return "RNA"
        if "C2'" in names:
            return "DNA"
    return None


def classify_structure(atoms):
    """Classify every chain and the structure as a whole.

    Returns
    -------
    chains : pandas.DataFrame
        One row per chain: residue counts per type and the chain ``label``
        (``DNA``, ``RNA``, ``DNA-RNA hybrid``, ``protein`` or ``other``).
    summary : str
        ``"DNA"``, ``"RNA"``, ``"both"`` or ``"neither"``.
    """
    rows = []
    # Group by chain ID rather than chain_iter(): HETATM records (e.g. waters) usually
    # come after all polymer chains, so chain_iter() would split one chain into two.
    for chain_id in pd.unique(atoms.chain_id):
        chain = atoms[atoms.chain_id == chain_id]
        counts = {"DNA": 0, "RNA": 0, "protein": 0, "other": 0}
        for res in struc.residue_iter(chain):
            counts[residue_type(res.res_name[0], res) or "other"] += 1

        n_nuc = counts["DNA"] + counts["RNA"]
        if counts["protein"] > n_nuc:
            label = "protein"
        elif counts["DNA"] and counts["RNA"]:
            label = "DNA-RNA hybrid"
        elif counts["DNA"]:
            label = "DNA"
        elif counts["RNA"]:
            label = "RNA"
        else:
            label = "other"

        rows.append({"chain": chain_id, **{f"n_{k}": v for k, v in counts.items()},
                     "label": label})

    chains = pd.DataFrame(rows)
    labels = set(chains["label"]) if len(chains) else set()
    has_dna = bool(labels & {"DNA", "DNA-RNA hybrid"})
    has_rna = bool(labels & {"RNA", "DNA-RNA hybrid"})
    summary = {(True, True): "both", (True, False): "DNA",
               (False, True): "RNA", (False, False): "neither"}[(has_dna, has_rna)]
    return chains, summary


def dna_chain_ids(chains):
    """Chain IDs whose label contains DNA (pure DNA or hybrid)."""
    return list(chains.loc[chains["label"].isin(["DNA", "DNA-RNA hybrid"]), "chain"])


def nucleic_chain_ids(chains):
    """Chain IDs labelled DNA, RNA or DNA-RNA hybrid."""
    return list(chains.loc[chains["label"].isin(["DNA", "RNA", "DNA-RNA hybrid"]), "chain"])
