"""R3 — chemical probing (SHAPE/DMS) consistency (external-validity tests).

The reactivity of each nucleotide (0..1) is compared against the paired /
unpaired state inferred geometrically from the structure: paired bases should
be unreactive, unpaired bases reactive. Mismatches beyond a fraction flag a
FAIL. The octagon geometry below gives every residue a spatial partner (i +/- 2),
so a fully consistent dataset passes and a fully contradictory one fails.
"""
import math

from rnvalidate.core.dataclasses import RnaInput, RnaPrediction, Residue
from rnvalidate.core.probing import check_probing

_R = 7.84  # radius chosen so adjacent C1' are ~6 A and i+/-2 are ~11.1 A (paired)


def _octagon():
    coords = {}
    for k in range(8):
        ang = 2 * math.pi * k / 8
        coords[k + 1] = (_R * math.cos(ang), _R * math.sin(ang), 0.0)
    residues = [Residue(resseq=i, resname="A", c1_coord=c) for i, c in coords.items()]
    return RnaPrediction(residues=residues)


def test_probing_ok_consistent():
    shape = {str(i): 0.10 for i in range(1, 9)}  # all paired -> all unreactive
    rna = RnaInput(pdb_id="t", structure=_octagon(), exp={"shape_reactivity": shape})
    assert check_probing(rna) is None


def test_probing_fail_contradictory():
    shape = {str(i): 0.90 for i in range(1, 9)}  # paired but reactive -> all mismatch
    rna = RnaInput(pdb_id="t", structure=_octagon(), exp={"shape_reactivity": shape})
    v = check_probing(rna)
    assert v is not None
    assert v.rule_id == "R3"
    assert v.severity == "FAIL"
    assert "probing_inconsistent" in v.reason


def test_probing_threshold_not_overfiring():
    # Only 1 of 8 positions contradicts -> fraction 0.125 < 0.40 -> PASS.
    shape = {str(i): 0.10 for i in range(1, 9)}
    shape["1"] = 0.90
    rna = RnaInput(pdb_id="t", structure=_octagon(), exp={"shape_reactivity": shape})
    assert check_probing(rna) is None


def test_probing_skip_without_data():
    rna = RnaInput(pdb_id="t", structure=_octagon(), exp={})
    v = check_probing(rna)
    assert v is not None
    assert v.rule_id == "R3"
    assert v.severity == "SKIP"
