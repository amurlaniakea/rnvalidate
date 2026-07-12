"""CLI end-to-end + external-validity + predictor-agnosticism (AC-8).

Invokes the tool under test via `python -m rnvalidate.cli` (no `pip install -e`
needed; the package is importable from the repo root). The only subprocess used
is the RNAValidate CLI itself -- no external predictor (AF3/RiboSphere/etc) is
ever launched from a test.
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
FIX = REPO / "rnvalidate" / "fixtures"


def _run(args):
    return subprocess.run(
        [sys.executable, "-m", "rnvalidate.cli", *args],
        capture_output=True, text=True, cwd=str(REPO),
    )


# --- E1: clean hairpin -> PASS, exact report keys -------------------------
def test_cli_e1_pass():
    r = _run(["check", "--in", str(FIX / "hairpin_clean.pdb"),
              "--exp", str(FIX / "fret_clean.json"), "--format", "json"])
    assert r.returncode == 0, r.stderr
    data = json.loads(r.stdout)
    assert isinstance(data, list) and len(data) == 1
    obj = data[0]
    assert obj["pdb_id"] == "hairpin_clean"
    assert obj["verdict"] == "PASS"
    assert {"rule_id", "severity", "reason"} <= set(obj["violations"][0].keys())
    # external-validity: no rule is FAIL for the known-good hairpin
    assert not any(v["severity"] == "FAIL" for v in obj["violations"])


# --- E2: high-ranked but hydrolyzing aptamer -> FAIL with F1 --------------
def test_cli_e2_fail():
    r = _run(["check", "--in", str(FIX / "aptamer_hydrolyzes.pdb"),
              "--exp", str(FIX / "fret_ok.json")])
    assert r.returncode == 1, r.stderr
    data = json.loads(r.stdout)
    obj = data[0]
    assert obj["verdict"] == "FAIL"
    ids = {v["rule_id"] for v in obj["violations"]}
    assert any(i.startswith("F") for i in ids)  # F1 present (tube false-positive)


# --- AC-9: markdown report to stdout --------------------------------------
def test_cli_markdown():
    r = _run(["check", "--in", str(FIX / "hairpin_clean.pdb"),
              "--exp", str(FIX / "fret_clean.json"), "--format", "md"])
    assert r.returncode == 0, r.stderr
    assert "RNAValidate" in r.stdout
    assert "|" in r.stdout


# --- dir mode + sidecar loading (covers _load_inputs dir branch) ----------
def test_cli_dir_with_sidecars():
    r = _run(["check", "--in", str(FIX), "--format", "json"])
    assert r.returncode == 1  # aptamer drags overall verdict to FAIL
    data = json.loads(r.stdout)
    assert len(data) == 2
    assert {"hairpin_clean", "aptamer_hydrolyzes"} <= {o["pdb_id"] for o in data}


# --- no --exp path (covers _load_exp(None)) -------------------------------
def test_cli_no_exp_does_not_crash():
    r = _run(["check", "--in", str(FIX / "hairpin_clean.pdb")])
    assert r.returncode in (0, 1), r.stderr
    json.loads(r.stdout)  # still valid JSON


# --- input error -> exit code 2 -------------------------------------------
def test_cli_bad_input_exit2():
    r = _run(["check", "--in", str(FIX / "does_not_exist.pdb")])
    assert r.returncode == 2


# --- unit coverage for the PDB parser + rules parse-error path -----------
def test_pdb_parse_error_propagates_to_r1_fail():
    from rnvalidate.core.pdb import parse_pdb
    from rnvalidate.core.rules import apply_rules
    from rnvalidate.core.dataclasses import RnaInput

    pred = parse_pdb("ATOM      1 C1'    G A   X   not-a-number here\n")
    assert pred.parse_error is not None
    v = apply_rules(RnaInput(pdb_id="x", structure=pred))
    assert v.verdict == "FAIL"
    assert any(x.rule_id == "R1" and x.severity == "FAIL" for x in v.violations)


def test_pdb_is_rna_detection():
    from rnvalidate.core.pdb import parse_pdb, is_rna

    rna_pred = parse_pdb("ATOM      1 C1'    A A   1       0.0 0.0 0.0\n")
    assert is_rna(rna_pred) is True
    prot_pred = parse_pdb("ATOM      1 CA     ALA A   1       0.0 0.0 0.0\n")
    assert is_rna(prot_pred) is False


# --- AC-8: no test ever launches an external predictor --------------------
def test_no_external_predictor_references():
    forbidden = [
        "alphafold", "af3", "rosettafold", "ribosphere", "colabfold",
        "esmfold", "openfold", "rfdiffusion", "af2",
    ]
    text = ""
    for f in (REPO / "tests").glob("test_*.py"):
        text += f.read_text(encoding="utf-8").lower()
    for word in forbidden:
        assert word not in text, f"predictor reference found in tests: {word}"
