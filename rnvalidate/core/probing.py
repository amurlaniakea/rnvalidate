"""R3 - Chemical probing (SHAPE / DMS).

La reactividad SHAPE/DMS reporta que nucleotidos estan apareados (baja
reactividad) o flexibles/no apareados (alta reactividad). Se compara con el
estado pareado/no-pareado PREDICHO por la estructura (heuristica de vecindad:
un residuo con un vecino espacial C1'-C1' cercano se considera apareado).

Si >SHAPE_MISMATCH_FRAC de las posiciones con dato contradicen la estructura
=> FAIL. Sin datos de probing => SKIP.

Contrato exp:
  "shape_reactivity": { "7": 0.85, "12": 0.10 }   (resseq -> reactividad 0..1)
"""
from __future__ import annotations

import math
from typing import Optional

from rnvalidate.core.constants import SHAPE_MISMATCH_FRAC, SHAPE_REACTIVE_THRESHOLD
from rnvalidate.core.dataclasses import RnaInput, RuleViolation

# Distancia C1'-C1' tipica de un par de bases Watson-Crick ~ 10.4 A.
_PAIR_DIST_A = 12.0


def _dist(a, b) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _is_paired(rna: RnaInput, resseq: int) -> Optional[bool]:
    """True si el residuo tiene un partner espacial (no secuencial) cercano."""
    coords = {r.resseq: r.c1_coord for r in rna.structure.residues
              if r.c1_coord is not None}
    if resseq not in coords:
        return None
    me = coords[resseq]
    for other, c in coords.items():
        if abs(other - resseq) <= 1:
            continue  # vecino secuencial, no cuenta como par
        if _dist(me, c) <= _PAIR_DIST_A:
            return True
    return False


def check_probing(rna: RnaInput) -> Optional[RuleViolation]:
    shape = rna.exp.get("shape_reactivity")
    if not shape:
        return RuleViolation("R3", "SKIP", "no_probing_data: sin SHAPE/DMS")
    if rna.structure.parse_error:
        return RuleViolation("R3", "SKIP", "pdb_unparseable")

    total = 0
    mismatches = 0
    for res_s, react in shape.items():
        try:
            resseq = int(res_s)
            react = float(react)
        except (ValueError, TypeError):
            continue
        paired = _is_paired(rna, resseq)
        if paired is None:
            continue
        total += 1
        reactive = react > SHAPE_REACTIVE_THRESHOLD
        # apareado deberia ser NO reactivo; no apareado deberia ser reactivo
        if paired and reactive:
            mismatches += 1
        elif (not paired) and (not reactive):
            mismatches += 1

    if total == 0:
        return RuleViolation("R3", "SKIP", "no_matching_positions")

    frac = mismatches / total
    if frac > SHAPE_MISMATCH_FRAC:
        return RuleViolation(
            "R3", "FAIL",
            f"probing_inconsistent: {mismatches}/{total} posiciones contradicen "
            f"la estructura (frac={frac:.2f} > {SHAPE_MISMATCH_FRAC})",
        )
    return None
