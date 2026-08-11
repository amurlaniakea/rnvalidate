# SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Constantes fisicas / umbrales para las reglas de RNAValidate.

NO estan tuneadas para que los fixtures pasen; se derivan de fisica/quimica
documentada:
- FRET_RMSD_TOL_A: FRET single-molecule tiene precision ~2-6 A en distancias
  30-80 A; un RMSD > 8 A entre distancias caidas de la estructura y las
  restricciones FRET indica inconsistencia estructural real, no ruido.
- EM_CORR_MIN: la correlacion cruzada mapa/modelo aceptable en cryo-EM de
  resolucion media es >0.5 (FSC/CC). Por debajo, el modelo no explica la densidad.
- SHAPE_MISMATCH_FRAC: SHAPE/DMS: nucleotidos apareados -> baja reactividad
  (<0.3), no apareados -> alta (>0.5). Si >40% de las posiciones con dato
  contradicen el estado pareado/no-pareado de la estructura, R3 FAIL.
- HYDROLYSIS_RISK_MAX: score derivado (0..1) de propension a auto-cleavage;
  >0.6 marca falso positivo de "validacion en tubo".
"""
from __future__ import annotations

FRET_RMSD_TOL_A: float = 8.0
EM_CORR_MIN: float = 0.5
SHAPE_MISMATCH_FRAC: float = 0.40
HYDROLYSIS_RISK_MAX: float = 0.60

# Reactividad SHAPE: umbral que separa "apareado" (bajo) de "no apareado" (alto).
SHAPE_REACTIVE_THRESHOLD: float = 0.4
