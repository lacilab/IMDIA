"""Loading PDB files into biotite AtomArrays."""

from pathlib import Path

import biotite.structure.io.pdb as pdb


def load_structure(path, model=1):
    """Load one model of a legacy ``.pdb`` file.

    Parameters
    ----------
    path : str or Path
        Path to the PDB file.
    model : int
        ``MODEL`` serial number to load, 1-based, as in the file. Biotite uses the
        same 1-based numbering, so no conversion is needed.

    Returns
    -------
    atoms : biotite.structure.AtomArray
        The selected model. For atoms with alternate locations, the altloc with
        the highest occupancy is kept.
    info : dict
        ``model_count``, ``model_used`` and a list of human-readable ``notes``.
    """
    path = Path(path)
    if path.suffix.lower() != ".pdb":
        raise ValueError(f"Only .pdb files are supported in this version, got {path.name}")

    pdb_file = pdb.PDBFile.read(str(path))
    model_count = pdb_file.get_model_count()
    if not 1 <= model <= model_count:
        raise ValueError(f"Model {model} requested, but the file has {model_count} model(s)")

    atoms = pdb_file.get_structure(model=model, altloc="occupancy")

    notes = []
    if model_count > 1:
        notes.append(
            f"File has {model_count} models; only model {model} was analysed, "
            f"the other {model_count - 1} were skipped."
        )
    return atoms, {"model_count": model_count, "model_used": model, "notes": notes}
