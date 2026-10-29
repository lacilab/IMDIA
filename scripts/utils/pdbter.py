"""
Use TER records to determine where strands end.
"""

from pathlib import Path


def ter_residues(pdb, model=1):
    """Residues immediately followed by a TER record in the given model.

    Returns
    -------
    set of tuple
        ``(chain_id, res_id, ins_code)`` of the last atom before each TER.
    """
    ter_after, current, last = set(), 1, None
    with open(Path(pdb)) as f:
        for line in f:
            if line.startswith("MODEL"):
                try:
                    current = int(line[10:14])
                except ValueError:
                    pass                     # malformed MODEL line: keep the current model
                last = None
            elif current != model:
                continue
            elif line.startswith(("ATOM", "HETATM")):
                try:
                    last = (line[21], int(line[22:26]), line[26].strip())
                except ValueError:
                    last = None
            elif line.startswith("TER") and last is not None:
                ter_after.add(last)
                last = None
    return ter_after


def residue_key(atoms, k):
    """Key of the residue that atom ``k`` belongs to."""
    return (atoms.chain_id[k], int(atoms.res_id[k]), atoms.ins_code[k])
