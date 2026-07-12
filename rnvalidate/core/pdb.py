"""Parser minimo de PDB para RNA (sin Biopython obligatorio).

Lee solo lineas ATOM/HETATM. Extrae residuos (nucleotidos) y la coordenada
del atomo C1' de cada residuo (usada por R1 FRET para calcular distancias).
El filtro SOLO LEE el PDB; no lo modifica ni invoca predictores.
"""
from __future__ import annotations

from typing import List

from rnvalidate.core.dataclasses import Residue, RnaPrediction

_RNA_RESNAMES = {"A", "U", "G", "C", "RA", "RU", "RG", "RC",
                 "ADE", "URA", "GUA", "CYT"}


def parse_pdb(text: str) -> RnaPrediction:
    """Parsea texto PDB a RnaPrediction. Marca parse_error si esta corrupto."""
    residues: List[Residue] = []
    n_atoms = 0
    seen = {}
    any_atom_line = False
    try:
        for line in text.splitlines():
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            any_atom_line = True
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            resseq = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            n_atoms += 1
            if resseq not in seen:
                r = Residue(resseq=resseq, resname=resname)
                seen[resseq] = r
                residues.append(r)
            if atom_name in ("C1'", "C1*"):
                seen[resseq].c1_coord = (x, y, z)
    except (ValueError, IndexError) as e:
        return RnaPrediction(residues=[], n_atoms=0,
                             parse_error=f"unparseable_pdb: {e}")

    if not any_atom_line:
        return RnaPrediction(residues=[], n_atoms=0,
                             parse_error="unparseable_pdb: no ATOM/HETATM records")

    residues.sort(key=lambda r: r.resseq)
    return RnaPrediction(residues=residues, n_atoms=n_atoms)


def is_rna(pred: RnaPrediction) -> bool:
    """True si al menos un residuo tiene un resname de nucleotido de RNA."""
    return any(r.resname.upper() in _RNA_RESNAMES for r in pred.residues)
