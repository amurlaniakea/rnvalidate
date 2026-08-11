# SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
# SPDX-License-Identifier: AGPL-3.0-or-later

"""R1 - Consistencia FRET.

Toma las distancias caidas de la estructura 3D (entre atomos C1' de pares de
residuos) y las compara con las restricciones FRET experimentales publicadas.
El RMSD de las distancias (estructura vs FRET) > tolerancia => FAIL.

Contrato exp:
  "fret_distances_angstrom": { "7-19": 38.2, "3-25": 51.0 }
donde la clave "i-j" son los resseq de los dos residuos etiquetados.

Sin datos FRET => SKIP (no FAIL). PDB no parseable => SKIP (R1_pdb lo cubre).
"""
from __future__ import annotations

import math
from typing import Dict, Optional

from rnvalidate.core.constants import FRET_RMSD_TOL_A
from rnvalidate.core.dataclasses import RnaInput, RuleViolation


def _dist(a, b) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _coord_by_resseq(rna: RnaInput) -> Dict[int, tuple]:
    return {r.resseq: r.c1_coord for r in rna.structure.residues
            if r.c1_coord is not None}


def check_fret(rna: RnaInput) -> Optional[RuleViolation]:
    fret = rna.exp.get("fret_distances_angstrom")
    if not fret:
        return RuleViolation("R1", "SKIP", "no_fret_data: sin restricciones FRET")
    if rna.structure.parse_error:
        return RuleViolation("R1", "SKIP", "pdb_unparseable: no se puede medir distancia")

    coords = _coord_by_resseq(rna)
    sq_errors = []
    for pair, exp_d in fret.items():
        try:
            i_s, j_s = pair.split("-")
            i, j = int(i_s), int(j_s)
        except ValueError:
            continue
        if i not in coords or j not in coords:
            continue
        struct_d = _dist(coords[i], coords[j])
        sq_errors.append((struct_d - float(exp_d)) ** 2)

    if not sq_errors:
        return RuleViolation("R1", "SKIP", "no_matching_pairs: pares FRET no en estructura")

    rmsd = math.sqrt(sum(sq_errors) / len(sq_errors))
    if rmsd > FRET_RMSD_TOL_A:
        return RuleViolation(
            "R1", "FAIL",
            f"fret_inconsistent: RMSD distancias={rmsd:.2f} A > tol {FRET_RMSD_TOL_A} A",
        )
    return None
