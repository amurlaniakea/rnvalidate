"""Render de veredictos a JSON, Markdown y SARIF 2.1.0."""

from __future__ import annotations

import json
from typing import List, Optional

from rnvalidate.core.dataclasses import Verdict, RnaInput, RuleViolation


def to_json(verdicts: List[Verdict]) -> str:
    out = []
    for v in verdicts:
        out.append({
            "pdb_id": v.pdb_id,
            "verdict": v.verdict,
            "violations": [
                {"rule_id": x.rule_id, "severity": x.severity, "reason": x.reason}
                for x in v.violations
            ],
        })
    return json.dumps(out, indent=2, ensure_ascii=False)


def to_markdown(verdicts: List[Verdict]) -> str:
    lines = ["# RNAValidate — reporte", "",
             "| PDB ID | Veredicto | Reglas |",
             "|--------|-----------|--------|"]
    for v in verdicts:
        if not v.violations:
            viol = "—"
        else:
            viol = "; ".join(f"{x.rule_id}({x.severity}): {x.reason}" for x in v.violations)
        lines.append(f"| `{v.pdb_id}` | **{v.verdict}** | {viol} |")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# SARIF 2.1.0 (feature 004)
# ---------------------------------------------------------------------------
# Reglas DURAS -> "error"; regla DERIVADA F1 -> "warning" (Opción B aprobada).
# El core mantiene F1 como severity=FAIL (MVP intacto); el mapeo es de
# presentación para GitHub Security, citando TRANSPARENCIA_RNAValidate.md.
LEVEL_BY_RULE = {
    "R1": "error", "R2": "error", "R3": "error", "R4": "error",
    "F1": "warning",
}

RULE_META = {
    "R1": {
        "name": "pdb_sanity",
        "shortDescription": {"text": "Estructura PDB parseable y enlaces coherentes"},
        "fullDescription": {"text": "Validación de sanity de la estructura PDB (enlaces, átomos)"},
        "helpUri": "https://github.com/amurlaniakea/rnvalidate",
    },
    "R2": {
        "name": "datos_exp",
        "shortDescription": {"text": "Consistencia con datos experimentales (FRET/cryo-EM/SHAPE)"},
        "fullDescription": {"text": "Contraste de la predicción con datos experimentales publicados"},
        "helpUri": "https://github.com/amurlaniakea/rnvalidate",
    },
    "R3": {
        "name": "plegamiento",
        "shortDescription": {"text": "Sin pseudoknots imposibles / topología coherente"},
        "fullDescription": {"text": "Sanity topológica del plegamiento"},
        "helpUri": "https://github.com/amurlaniakea/rnvalidate",
    },
    "R4": {
        "name": "energia",
        "shortDescription": {"text": "Energía en rango físico"},
        "fullDescription": {"text": "Energía de la estructura dentro de rango físico plausible"},
        "helpUri": "https://github.com/amurlaniakea/rnvalidate",
    },
    "F1": {
        "name": "hidrolisis",
        "shortDescription": {"text": "Propensión a hidrólisis cinética (DERIVADA, requiere half-life en tubo)"},
        "fullDescription": {"text": "Heurística de hidrólisis: señal de falso positivo en tubo, NO certificado"},
        "helpUri": "https://github.com/amurlaniakea/rnvalidate",
    },
}


def _driver_version() -> str:
    try:
        import importlib.metadata
        return importlib.metadata.version("rnvalidate")
    except Exception:
        return "0.1.0"


def _driver_rules() -> List[dict]:
    rules = []
    for rid in ("R1", "R2", "R3", "R4", "F1"):
        meta = RULE_META[rid]
        rules.append({
            "id": rid,
            "name": meta["name"],
            "shortDescription": meta["shortDescription"],
            "fullDescription": meta["fullDescription"],
            "helpUri": meta["helpUri"],
        })
    return rules


def _sarif_level(viol: "RuleViolation") -> str:
    """Nivel SARIF por regla (mapeo de presentación; core intacto)."""
    return LEVEL_BY_RULE.get(viol.rule_id, "error")


def to_sarif(verdicts: List[Verdict], inputs: Optional[List[RnaInput]] = None) -> str:
    """Renderiza veredictos como documento SARIF 2.1.0 (str JSON).

    Solo las violaciones con severity="FAIL" generan 'results'. Los skips
    (severity="SKIP") no producen result. El mapeo F1 FAIL -> "warning" es de
    presentación únicamente (ver LEVEL_BY_RULE). No se emite 'region'.
    """
    if inputs is None:
        inputs = []

    input_by_id: dict[str, RnaInput] = {m.pdb_id: m for m in inputs}

    results = []
    for v in verdicts:
        mol = input_by_id.get(v.pdb_id)
        uri = (mol and getattr(mol, "pdb_path", None)) or "--pdb"
        for viol in v.violations:
            if viol.severity != "FAIL":
                continue  # skips no generan result
            results.append({
                "ruleId": viol.rule_id,
                "level": _sarif_level(viol),
                "message": {"text": viol.reason},
                "locations": [
                    {
                        "physicalLocation": {
                            "artifactLocation": {
                                "uri": uri,
                                "properties": {
                                    "moleculeId": v.pdb_id,
                                    "ruleId": viol.rule_id,
                                },
                            }
                        }
                    }
                ],
            })

    doc = {
        "$schema": "https://json.schemastore.org/sarif-2.1.0.json",
        "version": "2.1.0",
        "runs": [
            {
                "tool": {
                    "driver": {
                        "name": "rnvalidate",
                        "version": _driver_version(),
                        "informationUri": "https://github.com/amurlaniakea/rnvalidate",
                        "rules": _driver_rules(),
                    }
                },
                "results": results,
            }
        ],
    }
    return json.dumps(doc, indent=2, ensure_ascii=False)
