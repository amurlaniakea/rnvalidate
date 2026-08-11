# SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
# SPDX-License-Identifier: AGPL-3.0-or-later

"""CLI end-to-end + external-validity + predictor-agnosticism (AC-8).

The CLI is exercised IN-PROCESS via Typer's CliRunner so that coverage is
measured for cli.py, report.py and the PDB parser. The only code under test is
RNAValidate itself -- no external structure predictor is ever launched from a
test (AC-8).
"""
import json
from pathlib import Path

from typer.testing import CliRunner

from rnvalidate.cli import app
from rnvalidate.core.pdb import parse_pdb, is_rna
from rnvalidate.core.rules import apply_rules
from rnvalidate.core.dataclasses import RnaInput, Residue, RnaPrediction

REPO = Path(__file__).resolve().parent.parent
FIX = REPO / "rnvalidate" / "fixtures"
runner = CliRunner()


def _pdb_line(atom, resname, resseq, x=0.0, y=0.0, z=0.0):
    """Build a fixed-column PDB ATOM line (resName at columns 18-20)."""
    return (
        "ATOM  "
        + "    1"
        + " "
        + f"{atom:>4}"
        + " "
        + f"{resname:>3}"
        + " A"
        + f"{resseq:>4}"
        + "    "
        + f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
    )


# --- E1: clean hairpin -> PASS, exact report keys -------------------------
def test_cli_e1_pass():
    r = runner.invoke(
        app,
        ["check", "--in", str(FIX / "hairpin_clean.pdb"),
         "--exp", str(FIX / "fret_clean.json"), "--format", "json"],
    )
    assert r.exit_code == 0, r.output
    data = json.loads(r.output)
    assert isinstance(data, list) and len(data) == 1
    obj = data[0]
    assert obj["pdb_id"] == "hairpin_clean"
    assert obj["verdict"] == "PASS"
    assert {"rule_id", "severity", "reason"} <= set(obj["violations"][0].keys())
    # external-validity: no rule is FAIL for the known-good hairpin
    assert not any(v["severity"] == "FAIL" for v in obj["violations"])


# --- E2: high-ranked but hydrolyzing aptamer -> FAIL with F1 --------------
def test_cli_e2_fail():
    r = runner.invoke(
        app,
        ["check", "--in", str(FIX / "aptamer_hydrolyzes.pdb"),
         "--exp", str(FIX / "fret_ok.json")],
    )
    assert r.exit_code == 1, r.output
    data = json.loads(r.output)
    obj = data[0]
    assert obj["verdict"] == "FAIL"
    ids = {v["rule_id"] for v in obj["violations"]}
    assert any(i.startswith("F") for i in ids)  # F1 present (tube false-positive)


# --- AC-9: markdown report to stdout --------------------------------------
def test_cli_markdown():
    r = runner.invoke(
        app,
        ["check", "--in", str(FIX / "hairpin_clean.pdb"),
         "--exp", str(FIX / "fret_clean.json"), "--format", "md"],
    )
    assert r.exit_code == 0, r.output
    assert "RNAValidate" in r.output
    assert "|" in r.output


# --- dir mode + sidecar loading (covers _load_inputs dir branch) ----------
def test_cli_dir_with_sidecars():
    r = runner.invoke(app, ["check", "--in", str(FIX), "--format", "json"])
    # aptamer_hydrolyzes has no matching *.json sidecar -> loads with no exp,
    # so it does NOT trip F1; both structures pass -> exit 0.
    assert r.exit_code == 0, r.output
    data = json.loads(r.output)
    assert len(data) == 2
    by = {o["pdb_id"]: o for o in data}
    assert "hairpin_clean" in by and "aptamer_hydrolyzes" in by
    assert by["hairpin_clean"]["verdict"] == "PASS"
    assert by["aptamer_hydrolyzes"]["verdict"] == "PASS"


# --- no --exp path (covers _load_exp(None)) -------------------------------
def test_cli_no_exp_does_not_crash():
    r = runner.invoke(app, ["check", "--in", str(FIX / "hairpin_clean.pdb")])
    assert r.exit_code in (0, 1), r.output
    json.loads(r.output)  # still valid JSON


# --- input error -> exit code 2 -------------------------------------------
def test_cli_bad_input_exit2():
    r = runner.invoke(app, ["check", "--in", str(FIX / "does_not_exist.pdb")])
    assert r.exit_code == 2


def test_cli_missing_in_option_exit2():
    r = runner.invoke(app, ["check"])
    assert r.exit_code == 2


# --- unit coverage for the PDB parser + rules parse-error path -----------
def test_pdb_parse_error_propagates_to_r1_fail():
    pred = parse_pdb("ATOM      1 C1'    G A   X   not-a-number here\n")
    assert pred.parse_error is not None
    v = apply_rules(RnaInput(pdb_id="x", structure=pred))
    assert v.verdict == "FAIL"
    assert any(x.rule_id == "R1" and x.severity == "FAIL" for x in v.violations)


def test_pdb_no_atoms_is_error():
    pred = parse_pdb("HEADER    RNA sample\n")
    assert pred.parse_error is not None
    assert pred.residues == []


def test_is_rna_detection():
    rna_pred = parse_pdb(_pdb_line("C1'", "A", 1))
    assert is_rna(rna_pred) is True
    prot_pred = parse_pdb(_pdb_line("CA", "ALA", 1))
    assert is_rna(prot_pred) is False


# --- AC-8: no test ever launches an external predictor --------------------
def test_no_external_predictor_references():
    # AC-8: tests must never launch an external structure predictor in a
    # subprocess. Scan every test source for known predictor tokens; the
    # definition line of this list is skipped below.
    forbidden = [
        "alphafold", "af3", "rosettafold", "ribosphere", "colabfold",
        "esmfold", "openfold", "rfdiffusion", "af2",
    ]
    for f in (REPO / "tests").glob("test_*.py"):
        for ln in f.read_text(encoding="utf-8").splitlines():
            low = ln.lower()
            # skip the list definition, its quoted token lines, and comments
            if "forbidden" in low or ln.strip().startswith('"') or ln.strip().startswith("#"):
                continue
            for word in forbidden:
                assert word not in low, f"predictor reference in {f.name}: {word}"
