"""R1 — FRET consistency (external-validity unit tests).

Each test builds an explicit RnaInput with independent ground truth
(coordinates + FRET restraints chosen by hand), never derived from the
checker's own output. A PASS means the structure-derived distance matches the
restraint within tolerance; a FAIL means the restraint is physically far off.
"""
from rnvalidate.core.dataclasses import RnaInput, RnaPrediction, Residue
from rnvalidate.core.fret import check_fret


def _struct(coords_by_resseq):
    residues = [
        Residue(resseq=i, resname="A", c1_coord=c)
        for i, c in coords_by_resseq.items()
    ]
    return RnaPrediction(residues=residues)


def test_fret_ok_distance_matches():
    rna = RnaInput(
        pdb_id="t",
        structure=_struct({1: (0.0, 0.0, 0.0), 2: (10.0, 0.0, 0.0)}),
        exp={"fret_distances_angstrom": {"1-2": 10.0}},
    )
    assert check_fret(rna) is None


def test_fret_fail_inconsistent():
    rna = RnaInput(
        pdb_id="t",
        structure=_struct({1: (0.0, 0.0, 0.0), 2: (10.0, 0.0, 0.0)}),
        exp={"fret_distances_angstrom": {"1-2": 40.0}},  # 30 A off -> impossible
    )
    v = check_fret(rna)
    assert v is not None
    assert v.rule_id == "R1"
    assert v.severity == "FAIL"
    assert "RMSD" in v.reason


def test_fret_skip_without_data():
    rna = RnaInput(pdb_id="t", structure=_struct({1: (0.0, 0.0, 0.0)}), exp={})
    v = check_fret(rna)
    assert v is not None
    assert v.rule_id == "R1"
    assert v.severity == "SKIP"


def test_fret_skip_unparseable_pdb():
    rna = RnaInput(
        pdb_id="t",
        structure=RnaPrediction(residues=[], parse_error="unparseable_pdb"),
        exp={"fret_distances_angstrom": {"1-2": 10.0}},
    )
    v = check_fret(rna)
    assert v is not None
    assert v.rule_id == "R1"
    assert v.severity == "SKIP"
    assert "unparseable" in v.reason
