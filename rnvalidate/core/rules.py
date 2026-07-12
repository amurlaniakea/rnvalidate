"""Registry de reglas y orquestador apply_rules().

Cada regla es una funcion pura (RnaInput) -> Optional[RuleViolation].
apply_rules recorre el registro, aplica las reglas, agrega violaciones/skips y
finaliza el veredicto (PASS/FAIL). PDB no parseable => R1 FAIL global y resto SKIP.
"""
from __future__ import annotations

from typing import Callable, List, Optional

from rnvalidate.core.dataclasses import RnaInput, RuleViolation, Verdict
from rnvalidate.core.designability import check_designability
from rnvalidate.core.em_corr import check_em_corr
from rnvalidate.core.fret import check_fret
from rnvalidate.core.hydrolysis import check_hydrolysis
from rnvalidate.core.probing import check_probing

# Orden estable: R1..R4, F1
RULES: List[Callable[[RnaInput], Optional[RuleViolation]]] = [
    check_fret,           # R1
    check_em_corr,        # R2
    check_probing,        # R3
    check_designability,  # R4
    check_hydrolysis,     # F1
]


def apply_rules(rna: RnaInput) -> Verdict:
    v = Verdict(pdb_id=rna.pdb_id, verdict="PASS")

    if rna.structure.parse_error:
        v.add(RuleViolation("R1", "FAIL", rna.structure.parse_error))
        v.finalize()
        return v

    for rule in RULES:
        v.add(rule(rna))
    v.finalize()
    return v
