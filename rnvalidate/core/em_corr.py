"""R2 - Correlacion cryo-EM.

Compara la correlacion cruzada mapa/modelo (cross-correlation / FSC) reportada
en los datos experimentales contra un minimo aceptable. Si no hay mapa de
densidad => SKIP. Correlacion baja => FAIL.

Contrato exp:
  "cryo_em_map": "/data/map.mrc"        (presencia => hay mapa)
  "cryo_em_correlation": 0.72           (CC modelo/mapa ya calculada, 0..1)

Este filtro NO carga el .mrc ni corre correlacion pesada (agnostico de GPU):
recibe la CC ya calculada por el usuario. Si hay mapa pero no CC => SKIP con
motivo (dato incompleto, no FAIL).
"""
from __future__ import annotations

import math
from typing import Optional

from rnvalidate.core.constants import EM_CORR_MIN
from rnvalidate.core.dataclasses import RnaInput, RuleViolation


def check_em_corr(rna: RnaInput) -> Optional[RuleViolation]:
    has_map = bool(rna.exp.get("cryo_em_map"))
    cc = rna.exp.get("cryo_em_correlation")

    if not has_map and cc is None:
        return RuleViolation("R2", "SKIP", "no_cryo_em_map: sin mapa de densidad")
    if cc is None:
        return RuleViolation("R2", "SKIP", "no_correlation_value: mapa presente sin CC")

    try:
        cc = float(cc)
    except (TypeError, ValueError):
        return RuleViolation("R2", "FAIL", "invalid_correlation: valor no numerico")
    if math.isnan(cc) or math.isinf(cc):
        return RuleViolation("R2", "FAIL", "invalid_correlation: NaN/inf")

    if cc < EM_CORR_MIN:
        return RuleViolation(
            "R2", "FAIL",
            f"low_em_correlation: CC={cc:.2f} < min {EM_CORR_MIN} "
            f"(modelo no explica la densidad)",
        )
    return None
