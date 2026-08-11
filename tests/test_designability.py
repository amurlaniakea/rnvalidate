# SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
# SPDX-License-Identifier: AGPL-3.0-or-later

"""R4 — designability / foldable topology (external-validity tests).

Heuristic, CPU-only: needs enough residues, no impossible backbone jump
between sequential C1' atoms (~5-9 A), and enough spatial pairing to be a real
folded motif rather than an extended strand.
"""
import math

from rnvalidate.core.dataclasses import RnaInput, RnaPrediction, Residue
from rnvalidate.core.designability import check_designability

_R = 7.84


def _octagon():
    residues = []
    for k in range(8):
        ang = 2 * math.pi * k / 8
        residues.append(
            Residue(resseq=k + 1, resname="A",
                    c1_coord=(_R * math.cos(ang), _R * math.sin(ang), 0.0))
        )
    return RnaPrediction(residues=residues)


def test_designability_ok_folded():
    rna = RnaInput(pdb_id="t", structure=_octagon())
    assert check_designability(rna) is None


def test_designability_fail_broken_backbone():
    residues = [
        Residue(resseq=1, resname="A", c1_coord=(0.0, 0.0, 0.0)),
        Residue(resseq=2, resname="A", c1_coord=(0.0, 20.0, 0.0)),  # 20 A jump
        Residue(resseq=3, resname="A", c1_coord=(0.0, 26.0, 0.0)),
        Residue(resseq=4, resname="A", c1_coord=(0.0, 32.0, 0.0)),
    ]
    rna = RnaInput(pdb_id="t", structure=RnaPrediction(residues=residues))
    v = check_designability(rna)
    assert v is not None
    assert v.rule_id == "R4"
    assert v.severity == "FAIL"
    assert "broken_backbone" in v.reason


def test_designability_fail_too_few_residues():
    residues = [
        Residue(resseq=1, resname="A", c1_coord=(0.0, 0.0, 0.0)),
        Residue(resseq=2, resname="A", c1_coord=(0.0, 6.0, 0.0)),
    ]
    rna = RnaInput(pdb_id="t", structure=RnaPrediction(residues=residues))
    v = check_designability(rna)
    assert v is not None
    assert v.rule_id == "R4"
    assert v.severity == "FAIL"
    assert "too_few_residues" in v.reason


def test_designability_skip_unparseable():
    pred = RnaPrediction(residues=[], parse_error="unparseable_pdb")
    rna = RnaInput(pdb_id="t", structure=pred)
    v = check_designability(rna)
    assert v is not None
    assert v.rule_id == "R4"
    assert v.severity == "SKIP"


def test_designability_fail_not_designable():
    # Extended strand: sequential C1' spaced 7 A (valid backbone) but no
    # spatial pairing (non-sequential residues > 12 A apart) -> not a folded
    # motif -> not_designable FAIL.
    residues = [
        Residue(resseq=i + 1, resname="A", c1_coord=(0.0, 7.0 * i, 0.0))
        for i in range(8)
    ]
    rna = RnaInput(pdb_id="t", structure=RnaPrediction(residues=residues))
    v = check_designability(rna)
    assert v is not None
    assert v.rule_id == "R4"
    assert v.severity == "FAIL"
    assert "not_designable" in v.reason
