"""Nucleic-acid classification and Watson-Crick base-pair detection (biotite-based).

Modules: io (PDB loading), classify (DNA/RNA/protein per chain), pairs (observed pairs),
infer (pairs inferred along duplex registers), diagram (text diagrams).
"""

from .io import load_structure
from .classify import classify_structure, nucleic_chain_ids
from .pairs import find_wc_pairs
from .diagram import pairs_diagram

__all__ = ["load_structure", "classify_structure", "nucleic_chain_ids",
           "find_wc_pairs", "pairs_diagram"]
