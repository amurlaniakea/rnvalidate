# SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
# SPDX-License-Identifier: AGPL-3.0-or-later

"""FUENTE UNICA de dataclasses para RNAValidate.

Todos los modulos importan desde aqui. NUNCA se redefine
RnaInput / RnaPrediction / Verdict / RuleViolation en otro modulo
(evita el bug de dataclass duplicado).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

Coord = Tuple[float, float, float]


@dataclass
class Residue:
    """Un residuo (nucleotido) de la estructura 3D."""
    resseq: int
    resname: str
    c1_coord: Optional[Coord] = None  # coordenada del atomo C1'


@dataclass
class RnaPrediction:
    """Estructura 3D RNA PREDICHA por un modelo externo (no por este filtro).

    Solo se LEE. Contiene los residuos parseados del PDB + conteo de atomos.
    """
    residues: List[Residue] = field(default_factory=list)
    n_atoms: int = 0
    parse_error: Optional[str] = None


@dataclass
class RnaInput:
    """Entrada del filtro: id + estructura 3D + datos experimentales + buffer."""
    pdb_id: str
    structure: RnaPrediction
    pdb_path: Optional[str] = None
    exp: Dict[str, Any] = field(default_factory=dict)
    buffer: Dict[str, Any] = field(default_factory=dict)


@dataclass
class RuleViolation:
    """Una regla violada (FAIL) o saltada por falta de dato (SKIP)."""
    rule_id: str          # "R1".."R4" | "F1"
    severity: str         # "FAIL" | "SKIP"
    reason: str


@dataclass
class Verdict:
    """Veredicto agregado para una estructura."""
    pdb_id: str
    verdict: str                       # "PASS" | "FAIL"
    violations: List[RuleViolation] = field(default_factory=list)

    def add(self, v: Optional[RuleViolation]) -> None:
        if v is not None:
            self.violations.append(v)

    def finalize(self) -> None:
        self.verdict = "FAIL" if any(v.severity == "FAIL" for v in self.violations) else "PASS"
