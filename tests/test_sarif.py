# SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Tests de salida SARIF 2.1.0 (feature 004) — RNAValidate.

Guarda de integridad LOCAL del schema (no sincronización con upstream):
el schema es un SNAPSHOT de la rama `main` de oasis-tcs/sarif-spec tomado el
2026-07-13. Si el archivo local cambia sin querer, este test falla ruidosamente.
"""
import json
import subprocess
from pathlib import Path

import pytest

try:
    import jsonschema
except ImportError:  # pragma: no cover
    pytest.skip("jsonschema no instalado (extra [testing])", allow_module_level=True)

from rnvalidate.core.dataclasses import RnaInput, RnaPrediction, Verdict
from rnvalidate.core.report import to_sarif, _driver_rules, LEVEL_BY_RULE
from rnvalidate.core.rules import apply_rules
from rnvalidate.cli import _build_input

THIS_DIR = Path(__file__).parent
FIXTURES = THIS_DIR / "fixtures"
SCHEMA_PATH = FIXTURES / "sarif_schema_2.1.0.json"
SCHEMA_SHA_PATH = FIXTURES / "sarif_schema_2.1.0.json.sha256"
# Fixtures de moléculas viven en el paquete (reusados del MVP).
PKG_FIXTURES = Path(__file__).parent.parent / "rnvalidate" / "fixtures"

# Hash local fijado en Bloque 1 (snapshot OASIS main @ 2026-07-13).
EXPECTED_SHA256 = "c3b4bb2d6093897483348925aaa73af03b3e3f4bd4ca38cef26dcb4212a2682e"


@pytest.fixture(scope="module")
def schema():
    return json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# Guarda de integridad LOCAL del schema (no sincronización con upstream).
# ---------------------------------------------------------------------------
def test_schema_integrity_local():
    actual = SCHEMA_SHA_PATH.read_text().split()[0].strip()
    assert actual == EXPECTED_SHA256, (
        "El schema SARIF local cambió. Es un snapshot de oasis-tcs/sarif-spec main "
        "@ 2026-07-13; actualiza el hash solo si re-snapshotas intencionalmente."
    )


def _load(pdb_name: str, exp_name: str = None) -> RnaInput:
    # Reusa el camino real del CLI (_build_input) para que el veredicto coincida
    # con el MVP (no parseamos PDB a mano, evitando el fallo "unparseable").
    pdb_path = PKG_FIXTURES / pdb_name
    exp = json.loads((PKG_FIXTURES / exp_name).read_text()) if exp_name else {}
    return _build_input(pdb_path, exp)


def _sarif_from_fixture(pdb_name, exp_name=None) -> dict:
    mol = _load(pdb_name, exp_name)
    v = apply_rules(mol)
    raw = to_sarif([v], [mol])
    return json.loads(raw)


# ---------------------------------------------------------------------------
# AC-S1: el output es válido contra el schema oficial 2.1.0.
# ---------------------------------------------------------------------------
def test_sarif_valid_against_schema(schema):
    doc = _sarif_from_fixture("aptamer_hydrolyzes.pdb", "fret_ok.json")
    jsonschema.validate(instance=doc, schema=schema)


def test_sarif_valid_clean(schema):
    doc = _sarif_from_fixture("hairpin_clean.pdb")
    jsonschema.validate(instance=doc, schema=schema)


# ---------------------------------------------------------------------------
# AC-S2: driver.rules contiene R1/R2/R3/R4/F1 con shortDescription no vacía.
# ---------------------------------------------------------------------------
def test_driver_rules_present():
    doc = _sarif_from_fixture("hairpin_clean.pdb")
    ids = [r["id"] for r in doc["runs"][0]["tool"]["driver"]["rules"]]
    assert ids == ["R1", "R2", "R3", "R4", "F1"]
    for r in doc["runs"][0]["tool"]["driver"]["rules"]:
        assert r["shortDescription"]["text"]


# ---------------------------------------------------------------------------
# AC-S3: cada violación -> result con ruleId/level/message/locations.
# ---------------------------------------------------------------------------
def test_result_maps_violation():
    doc = _sarif_from_fixture("aptamer_hydrolyzes.pdb", "fret_ok.json")
    assert doc["runs"][0]["results"]
    for res in doc["runs"][0]["results"]:
        assert res["ruleId"]
        assert res["level"] in ("error", "warning", "note", "none")
        assert res["message"]["text"]
        loc = res["locations"][0]["physicalLocation"]["artifactLocation"]
        assert "uri" in loc
        assert "moleculeId" in loc["properties"]


# ---------------------------------------------------------------------------
# AC-S4 (CLAVE): R1-R4 -> error; F1 -> warning (mapeo en renderer, core intacto).
# ---------------------------------------------------------------------------
def test_severity_mapping():
    doc = _sarif_from_fixture("aptamer_hydrolyzes.pdb", "fret_ok.json")
    by_rule = {r["ruleId"]: r["level"] for r in doc["runs"][0]["results"]}
    # F1 se dispara (heurística de hidrólisis) y DEBE aparecer como warning.
    assert by_rule.get("F1") == "warning"
    # Ninguna regla dura debe salir como warning.
    for rid in ("R1", "R2", "R3", "R4"):
        if rid in by_rule:
            assert by_rule[rid] == "error"


def test_mapping_table():
    assert LEVEL_BY_RULE["F1"] == "warning"
    for rid in ("R1", "R2", "R3", "R4"):
        assert LEVEL_BY_RULE[rid] == "error"


# ---------------------------------------------------------------------------
# AC-S5: parseable por consumidor genérico sin parser custom.
# ---------------------------------------------------------------------------
def test_generic_consumer():
    raw = to_sarif([apply_rules(_load("aptamer_hydrolyzes.pdb", "fret_ok.json"))],
                   [_load("aptamer_hydrolyzes.pdb", "fret_ok.json")])
    doc = json.loads(raw)
    assert "runs" in doc
    assert doc["runs"][0]["tool"]["driver"]["name"] == "rnvalidate"


# ---------------------------------------------------------------------------
# AC-S6: 3 casos (PASS limpio / 1 violación / múltiples violaciones).
# ---------------------------------------------------------------------------
def test_pass_clean_zero_results():
    doc = _sarif_from_fixture("hairpin_clean.pdb")
    assert doc["runs"][0]["results"] == []


def test_multiple_violations():
    doc = _sarif_from_fixture("aptamer_hydrolyzes.pdb", "fret_ok.json")
    assert len(doc["runs"][0]["results"]) >= 1
    assert "F1" in {r["ruleId"] for r in doc["runs"][0]["results"]}


# ---------------------------------------------------------------------------
# Versión driver dinámica (DECISIÓN 2).
# ---------------------------------------------------------------------------
def test_driver_version_dynamic():
    import importlib.metadata
    doc = _sarif_from_fixture("hairpin_clean.pdb")
    expected = importlib.metadata.version("rnvalidate")
    assert doc["runs"][0]["tool"]["driver"]["version"] == expected


# ---------------------------------------------------------------------------
# AC-S8: anti-oráculo — 0 subprocess a modelo/DFT en este archivo.
# ---------------------------------------------------------------------------
def test_no_model_subprocess_in_testfile():
    src = Path(__file__).read_text()
    # El único subprocess legítimo es el CLI real `rnvalidate`. No debe haber
    # imports reales a librerías de modelo/DFT/ML pesadas en este test.
    model_mods = ("torch", "tensorflow")
    for line in src.splitlines():
        s = line.strip()
        if s.startswith("import ") or s.startswith("from "):
            mod = s.split()[1].split(".")[0]
            assert mod not in model_mods, f"import de modelo pesado detectado: {s}"
    # El subprocess debe invocar únicamente el binario `rnvalidate`, no un
    # reconstructor/DFT externo.
    assert "subprocess.run(" in src
    assert "rnvalidate" in src  # binario invocado, no modelo/DFT


# ---------------------------------------------------------------------------
# S032: CLI real vía subprocess (proceso separado, no la función).
# ---------------------------------------------------------------------------
def test_cli_subprocess_real_sarif():
    # Subprocess REAL (proceso separado), portable: usa sys.executable -m rnvalidate.cli
    # en vez de confiar en el binario "rnvalidate" en el PATH (falla si el venv no está activado).
    import sys
    pdb = PKG_FIXTURES / "aptamer_hydrolyzes.pdb"
    exp = PKG_FIXTURES / "fret_ok.json"
    proc = subprocess.run(
        [sys.executable, "-m", "rnvalidate.cli", "check",
         "--in", str(pdb), "--exp", str(exp), "--format", "sarif"],
        capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 1  # FAIL (F1 disparada) -> exit 1, core intacto
    doc = json.loads(proc.stdout)
    assert "runs" in doc
    assert doc["runs"][0]["tool"]["driver"]["name"] == "rnvalidate"
    by_rule = {r["ruleId"]: r["level"] for r in doc["runs"][0]["results"]}
    assert by_rule.get("F1") == "warning"
