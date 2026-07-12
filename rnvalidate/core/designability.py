"""R4 - Designabilidad.

Comprueba que la estructura objetivo sea foldable/designable como heuristica
de topologia (no folding real, CPU-only):
- Debe haber suficientes residuos (una horquilla real tiene >= 4 nt).
- La cadena no debe presentar saltos geometricos imposibles entre residuos
  consecutivos (distancia C1'-C1' de vecinos secuenciales fuera de rango fisico
  ~ 5-9 A indica una cadena rota / no plegable).
- La densidad de contactos (fraccion de residuos con al menos un par espacial)
  debe superar un minimo: una estructura sin ningun apareamiento no es una
  estructura RNA designable, es una hebra extendida.

No designable => FAIL. PDB no parseable => SKIP (lo cubre R1_pdb en rules).
"""
from __future__ import annotations

import math
from typing import Optional

from rnvalidate.core.dataclasses import RnaInput, RuleViolation

_MIN_RESIDUES = 4
_BACKBONE_MIN_A = 4.0
_BACKBONE_MAX_A = 9.0
_PAIR_DIST_A = 12.0
_MIN_PAIRED_FRAC = 0.15


def _dist(a, b) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def check_designability(rna: RnaInput) -> Optional[RuleViolation]:
    if rna.structure.parse_error:
        return RuleViolation("R4", "SKIP", "pdb_unparseable")

    res = [r for r in rna.structure.residues if r.c1_coord is not None]
    if len(res) < _MIN_RESIDUES:
        return RuleViolation(
            "R4", "FAIL",
            f"too_few_residues: {len(res)} < {_MIN_RESIDUES} (no es estructura plegable)",
        )

    # cadena rota: salto backbone imposible entre vecinos secuenciales
    for a, b in zip(res, res[1:]):
        if b.resseq - a.resseq == 1:
            d = _dist(a.c1_coord, b.c1_coord)
            if d < _BACKBONE_MIN_A or d > _BACKBONE_MAX_A:
                return RuleViolation(
                    "R4", "FAIL",
                    f"broken_backbone: C1'-C1' {a.resseq}->{b.resseq} = {d:.2f} A "
                    f"fuera de [{_BACKBONE_MIN_A}, {_BACKBONE_MAX_A}] (no plegable)",
                )

    # densidad de contactos: al menos algun apareamiento espacial
    paired = 0
    for i, r in enumerate(res):
        for j, o in enumerate(res):
            if i == j or abs(o.resseq - r.resseq) <= 1:
                continue
            if _dist(r.c1_coord, o.c1_coord) <= _PAIR_DIST_A:
                paired += 1
                break
    frac = paired / len(res)
    if frac < _MIN_PAIRED_FRAC:
        return RuleViolation(
            "R4", "FAIL",
            f"not_designable: fraccion apareada={frac:.2f} < {_MIN_PAIRED_FRAC} "
            f"(hebra extendida, sin motivo estable)",
        )
    return None
