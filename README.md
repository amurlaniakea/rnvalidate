# RNAValidate

> Perimeter validator / ranker for AI-predicted 3D RNA structures against experimental data.

RNAValidate is a CPU-only, predictor-agnostic "linter" for 3D RNA structures. You
feed it a structure predicted by **your** model (AlphaFold-3, RiboSphere,
RoseTTAFold-RNA, …) plus experimental data (FRET, cryo-EM, SHAPE/DMS), and it
tells you, with a per-rule reason, whether that structure is consistent with the
data and — critically — whether it will survive in the tube or degrade before
you can measure it (a "tube false-positive").

It does **not** predict or design structures, and it never invokes any external
predictor. It only *reads* the PDB/JSON you produce.

## Why

Predictors and designers are plentiful; an open-source, reproducible layer that
**audits** a predicted 3D RNA structure against experimental evidence is not.
Before spending thousands of USD synthesizing an aptamer that ranks high
in-silico but hydrolyzes in the tube, RNAValidate filters it.

## Rules

| ID | Check | FAIL condition |
|----|-------|----------------|
| R1 | FRET consistency | RMSD of structure-derived vs FRET distances > 8 Å |
| R2 | cryo-EM correlation | map/model CC < 0.5 (skipped without a map) |
| R3 | Chemical probing (SHAPE/DMS) | >40% of probed positions contradict the structure |
| R4 | Designability | too few residues / broken backbone / no foldable motif |
| F1 | Hydrolysis propensity (DERIVED) | score > 0.60 → likely tube false-positive |

R2/R3 are optional: without the corresponding data they are *skipped*, never
*FAIL*. Thresholds are derived from documented physics/chemistry, not tuned to
the fixtures.

> **Transparency.** F1 is a DERIVED rule (risk R5), proposed by the author and
> not extracted from a single paper. It needs experimental validation against
> measured degradation before promotion to a hard rule.

## Install

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install typer numpy pytest pytest-cov
# or, from this repo:
pip install -e .
```

## Usage

```bash
rnvalidate check --in structure.pdb --exp exp.json --format json
rnvalidate check --in fixtures/ --format markdown   # directory + sidecars
```

Exit codes: `0` = all PASS, `1` = at least one FAIL (tube false-positive /
inconsistency), `2` = input error.

## Tests

```bash
pytest --cov=rnvalidate
```

The suite is CPU-only and fixture-driven. Fixtures use real PDB structures plus
published/example experimental data (external-validity, not circular). No test
launches an external predictor (AC-8).

## License

AGPL-3.0-or-later © 2026 Pedro Sordo Martínez — amurlaniakea@gmail.com.
See [LICENSE](LICENSE).
