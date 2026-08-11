# SPDX-FileCopyrightText: 2026 Pedro Sordo Martínez <amurlaniakea@gmail.com>
# SPDX-License-Identifier: AGPL-3.0-or-later

"""CLI de RNAValidate.

Uso:
  rnvalidate check --in <archivo.pdb | directorio> [--exp <exp.json>] [--format json|md]

El input es PDB + JSON producidos por el usuario; el filtro SOLO LEE.
Las reglas son puras y NO invocan ningun predictor externo (agnostico de modelo).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import List, Optional

import typer

from rnvalidate.core.dataclasses import RnaInput
from rnvalidate.core.pdb import parse_pdb
from rnvalidate.core.report import to_json, to_markdown, to_sarif
from rnvalidate.core.rules import apply_rules

app = typer.Typer(help="Validador/ranking perimetral de estructuras 3D RNA vs datos experimentales")


@app.callback()
def _root() -> None:
    """Validador/ranking perimetral de estructuras 3D RNA vs datos experimentales."""


def _build_input(pdb_path: Path, exp_data: dict) -> RnaInput:
    text = pdb_path.read_text(encoding="utf-8")
    structure = parse_pdb(text)
    return RnaInput(
        pdb_id=exp_data.get("pdb_id", pdb_path.stem),
        structure=structure,
        pdb_path=str(pdb_path),
        exp=exp_data.get("exp", {}),
        buffer=exp_data.get("buffer", {}),
    )


def _load_exp(exp_path: Optional[Path]) -> dict:
    if exp_path is None:
        return {}
    return json.loads(exp_path.read_text(encoding="utf-8"))


def _load_inputs(in_path: Optional[Path], exp_path: Optional[Path]) -> List[RnaInput]:
    if not in_path:
        raise typer.BadParameter("Usa --in <archivo.pdb|dir>")
    if in_path.is_dir():
        inputs = []
        for p in sorted(in_path.glob("*.pdb")):
            sidecar = p.with_suffix(".json")
            exp_data = _load_exp(sidecar) if sidecar.exists() else {}
            inputs.append(_build_input(p, exp_data))
        return inputs
    exp_data = _load_exp(exp_path)
    return [_build_input(in_path, exp_data)]


@app.command()
def check(
    in_path: Path = typer.Option(None, "--in", help="archivo PDB o directorio"),
    exp: Path = typer.Option(None, "--exp", help="JSON de datos experimentales"),
    format: str = typer.Option("json", "--format", help="json | md | sarif"),
):
    """Valida una estructura 3D RNA contra datos experimentales y reglas duras."""
    try:
        rnas = _load_inputs(in_path, exp)
    except (json.JSONDecodeError, KeyError, OSError, typer.BadParameter) as e:
        typer.echo(f"ERROR de entrada: {e}", err=True)
        raise typer.Exit(code=2)

    verdicts = [apply_rules(r) for r in rnas]
    if format == "sarif":
        out = to_sarif(verdicts, rnas)
    elif format == "md":
        out = to_markdown(verdicts)
    else:
        out = to_json(verdicts)
    typer.echo(out)

    if any(v.verdict == "FAIL" for v in verdicts):
        raise typer.Exit(code=1)
    raise typer.Exit(code=0)


def main():
    app()


if __name__ == "__main__":
    main()
