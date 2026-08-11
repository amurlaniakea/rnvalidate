# SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
# SPDX-License-Identifier: AGPL-3.0-or-later

"""F1 - Propension de hidrolisis cinetica / falsos positivos en tubo.

Feature DERIVADA (regla propuesta, marcada como riesgo R5): estima la propension
del diseno a hidrolisis espontanea (auto-cleavage del backbone fosfodiester)
ANTES de poder medirse, y marca los "falsos positivos de validacion en tubo":
estructuras que rankean alto in-silico pero se degradan en el tubo.

Heuristica real y ejecutable (CPU-only), basada en quimica del RNA documentada.
El ataque en linea del 2'-OH sobre el fosfato adyacente (transesterificacion) es
la ruta dominante de hidrolisis de RNA y depende de:
  1. Residuos NO apareados / en bucles flexibles  -> el 2'-OH puede alinearse
     en-line (geometria SN2). Los tramos pareados (helice A) estan protegidos.
  2. Composicion pirimidinica (U, C) -> los enlaces UpA / CpA son puntos
     calientes de hidrolisis conocidos.
  3. pH del buffer > 7 -> desprotonacion del 2'-OH acelera el ataque
     (dependencia log-lineal de la constante de velocidad con [OH-]).
  4. Cationes divalentes (Mg2+) -> catalizan la transesterificacion en linea.

Score derivado (0..1):
  risk = w_loop*loop_frac + w_pyr*pyrimidine_loop_frac + w_pH*pH_factor
         + w_mg*mg_factor
score > HYDROLYSIS_RISK_MAX => F1 FAIL (falso positivo de tubo).
"""
from __future__ import annotations

import math
from typing import Optional

from rnvalidate.core.constants import HYDROLYSIS_RISK_MAX
from rnvalidate.core.dataclasses import RnaInput, RuleViolation

_PAIR_DIST_A = 12.0
_W_LOOP = 0.40
_W_PYR = 0.25
_W_PH = 0.20
_W_MG = 0.15
_PYRIMIDINES = {"U", "C", "RU", "RC", "URA", "CYT"}


def _dist(a, b) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _unpaired_residues(rna: RnaInput):
    res = [r for r in rna.structure.residues if r.c1_coord is not None]
    unpaired = []
    for r in res:
        paired = False
        for o in res:
            if o.resseq == r.resseq or abs(o.resseq - r.resseq) <= 1:
                continue
            if _dist(r.c1_coord, o.c1_coord) <= _PAIR_DIST_A:
                paired = True
                break
        if not paired:
            unpaired.append(r)
    return res, unpaired


def hydrolysis_risk(rna: RnaInput) -> float:
    """Score de propension a hidrolisis, 0..1. Expuesto para tests unitarios."""
    res, unpaired = _unpaired_residues(rna)
    if not res:
        return 0.0

    loop_frac = len(unpaired) / len(res)
    pyr_loop = sum(1 for r in unpaired if r.resname.upper() in _PYRIMIDINES)
    pyrimidine_loop_frac = pyr_loop / len(res)

    pH = float(rna.buffer.get("pH", 7.0))
    # pH 7 -> 0 ; pH 9 -> 1 (saturado). Por debajo de 7, sin contribucion.
    pH_factor = max(0.0, min(1.0, (pH - 7.0) / 2.0))

    mg_mM = float(rna.buffer.get("mg_mM", 0.0))
    # 0 mM -> 0 ; >=10 mM -> 1.
    mg_factor = max(0.0, min(1.0, mg_mM / 10.0))

    risk = (_W_LOOP * loop_frac + _W_PYR * pyrimidine_loop_frac
            + _W_PH * pH_factor + _W_MG * mg_factor)
    return round(min(1.0, risk), 3)


def check_hydrolysis(rna: RnaInput) -> Optional[RuleViolation]:
    if rna.structure.parse_error:
        return RuleViolation("F1", "SKIP", "pdb_unparseable")
    risk = hydrolysis_risk(rna)
    if risk > HYDROLYSIS_RISK_MAX:
        return RuleViolation(
            "F1", "FAIL",
            f"tube_false_positive: propension_hidrolisis={risk} > {HYDROLYSIS_RISK_MAX} "
            f"(regla propuesta R5, requiere validacion experimental)",
        )
    return None
