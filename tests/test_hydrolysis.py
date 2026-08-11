# SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
# SPDX-License-Identifier: AGPL-3.0-or-later

"""F1 — hydrolysis propensity / tube false-positive (external-validity tests).

DERIVED rule (risk R5): estimates spontaneous backbone phosphodiester
hydrolysis from loop fraction, pyrimidine loops, buffer pH and Mg2+. A stable,
well-paired design at neutral pH scores ~0; an exposed, single-stranded design
at alkaline pH with Mg2+ scores high and FAILs (tube false-positive).
"""
import math

from rnvalidate.core.dataclasses import RnaInput, RnaPrediction, Residue
from rnvalidate.core.hydrolysis import check_hydrolysis, hydrolysis_risk


def _octagon():
    residues = []
    R = 7.84
    for k in range(8):
        ang = 2 * math.pi * k / 8
        residues.append(
            Residue(resseq=k + 1, resname="A",
                    c1_coord=(R * math.cos(ang), R * math.sin(ang), 0.0))
        )
    return RnaPrediction(residues=residues)


def _extended_unpaired(n=10):
    residues = [
        Residue(resseq=i + 1, resname="U", c1_coord=(0.0, 30.0 * i, 0.0))
        for i in range(n)
    ]
    return RnaPrediction(residues=residues)


def test_hydrolysis_ok_stable():
    rna = RnaInput(pdb_id="t", structure=_octagon(),
                   buffer={"pH": 7.0, "mg_mM": 0.0})
    assert hydrolysis_risk(rna) == 0.0
    assert check_hydrolysis(rna) is None


def test_hydrolysis_fail_exposed_alkaline():
    rna = RnaInput(pdb_id="t", structure=_extended_unpaired(),
                   buffer={"pH": 9.0, "mg_mM": 10.0})
    risk = hydrolysis_risk(rna)
    assert 0.0 <= risk <= 1.0
    assert risk > 0.6  # documents the threshold boundary
    v = check_hydrolysis(rna)
    assert v is not None
    assert v.rule_id == "F1"
    assert v.severity == "FAIL"
    assert "tube_false_positive" in v.reason


def test_hydrolysis_skip_unparseable():
    pred = RnaPrediction(residues=[], parse_error="unparseable_pdb")
    rna = RnaInput(pdb_id="t", structure=pred, buffer={"pH": 9.0, "mg_mM": 10.0})
    v = check_hydrolysis(rna)
    assert v is not None
    assert v.rule_id == "F1"
    assert v.severity == "SKIP"
