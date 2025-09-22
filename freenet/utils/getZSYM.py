"""
freenet.utils.getZSYM
---------------------

Parse a PDB or GRO file and return (or write) its element symbols
in free‐format ZSYM form (one symbol per line).
"""

from pathlib import Path
from typing import List, Optional
import MDAnalysis as mda

__all__ = ["getZSYM"]

def getZSYM(
    input_path: str,
    output_path: Optional[str] = None,
) -> List[str]:
    """
    Read a .pdb or .gro file, extract the element symbols,
    and either return them as a list or write them to output_path.

    Parameters
    ----------
    input_path : str
        Path to your .pdb or .gro file
    output_path : Optional[str]
        If given, write symbols (one per line) to this file.

    Returns
    -------
    List[str]
        The list of element symbols (capitalized).
    """
    p = Path(input_path)

    syms = _parse_structure(p)

    if output_path:
        with open(output_path, "w") as f:
            f.write("\n".join(syms))
    return syms

def _parse_structure(structure_path: Path) -> List[str]:

    symbols: List[str] = []

    u = mda.Universe(structure_path)

    protein = u.select_atoms("protein")

    for atom in protein:
    # atom.type is topology‐dependent; if missing, you might fall back to atom.element
        atype = atom.type if atom.type is not None else atom.element

        symbols.append(atype)


    return symbols