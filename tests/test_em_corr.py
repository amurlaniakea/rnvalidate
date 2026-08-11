# SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
# SPDX-License-Identifier: AGPL-3.0-or-later

"""R2 — cryo-EM cross-correlation (external-validity unit tests).

The checker reads a pre-computed map/model correlation (0..1) supplied by the
user; it never loads a .mrc or runs GPU correlation (predictor-agnostic).
"""
import math

from rnvalidate.core.dataclasses import RnaInput, RnaPrediction
from rnvalidate.core.em_corr import check_em_corr


def _empty():
    return RnaPrediction(residues=[])


def test_em_ok_good_correlation():
    rna = RnaInput(
        pdb_id="t",
        structure=_empty(),
        exp={"cryo_em_map": "/data/map.mrc", "cryo_em_correlation": 0.82},
    )
    assert check_em_corr(rna) is None


def test_em_fail_low_correlation():
    rna = RnaInput(
        pdb_id="t",
        structure=_empty(),
        exp={"cryo_em_map": "/data/map.mrc", "cryo_em_correlation": 0.20},
    )
    v = check_em_corr(rna)
    assert v is not None
    assert v.rule_id == "R2"
    assert v.severity == "FAIL"
    assert "low_em_correlation" in v.reason


def test_em_skip_without_map():
    rna = RnaInput(pdb_id="t", structure=_empty(), exp={})
    v = check_em_corr(rna)
    assert v is not None
    assert v.rule_id == "R2"
    assert v.severity == "SKIP"


def test_em_skip_map_without_correlation():
    rna = RnaInput(pdb_id="t", structure=_empty(), exp={"cryo_em_map": "/data/map.mrc"})
    v = check_em_corr(rna)
    assert v is not None
    assert v.rule_id == "R2"
    assert v.severity == "SKIP"


def test_em_fail_nan_correlation():
    rna = RnaInput(
        pdb_id="t",
        structure=_empty(),
        exp={"cryo_em_map": "/data/map.mrc", "cryo_em_correlation": float("nan")},
    )
    v = check_em_corr(rna)
    assert v is not None
    assert v.severity == "FAIL"
    assert "NaN" in v.reason


def test_em_fail_invalid_correlation():
    rna = RnaInput(
        pdb_id="t",
        structure=_empty(),
        exp={"cryo_em_map": "/data/map.mrc", "cryo_em_correlation": "abc"},
    )
    v = check_em_corr(rna)
    assert v is not None
    assert v.severity == "FAIL"
    assert "invalid_correlation" in v.reason
