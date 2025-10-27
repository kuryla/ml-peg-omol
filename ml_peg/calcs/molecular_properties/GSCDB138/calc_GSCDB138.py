"""Integrate GSCDB138 benchmark with per‑subset results."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import mlipx
from mlipx.abc import NodeWithCalculator
import pytest
import zntrack

from ml_peg.calcs.utils.utils import chdir
from ml_peg.models.get_models import load_models
from ml_peg.models.models import current_models


MODELS = load_models(current_models)


OUT_PATH = Path(__file__).parent / "outputs"


def get_external_script() -> Optional[Path]:
    """Locate external GSCDB138 evaluator if available."""
    env = os.environ.get("GSCDB138_SCRIPT")
    if env:
        p = Path(env)
        return p if p.exists() else None
    p = Path("/Users/Lilyes/Documents/GitHub/benchmarks-mp/gscdb138/eval_gscdb138_field.py")
    return p if p.exists() else None


class GSCDB138Benchmark(zntrack.Node):
    """
    Bootstrap node for GSCDB138.

    Standardises outputs under outputs/<model>/ with either a single CSV
    'gscdb138_results.csv' containing columns: subset,entry_id,ref,pred
    or multiple per‑subset CSVs named <subset>.csv with the same columns.

    This node doesn't run the heavy evaluation itself; instead, it prepares
    the output directory and prints clear instructions and expected formats.
    """

    model: NodeWithCalculator = zntrack.deps()
    model_name: str = zntrack.params()

    def run(self):
        model_out = OUT_PATH / self.model_name
        model_out.mkdir(parents=True, exist_ok=True)

        script = get_external_script()
        lines = [
            f"[GSCDB138] Prepared for model: {self.model_name}",
            f"Write results to: {model_out}",
            "Accepted:",
            " - Single: gscdb138_results.csv (subset,entry_id,ref,pred)",
            " - Multiple: <subset>.csv (subset,entry_id,ref,pred)",
        ]
        if script is not None:
            lines.append(f"External script (optional): {script}")
        else:
            lines.append("Set GSCDB138_SCRIPT to your evaluation script path.")
        print("\n".join(lines))


def build_project(repro: bool = False) -> None:
    project = mlipx.Project()
    nodes: dict[str, GSCDB138Benchmark] = {}
    for model_name, model in MODELS.items():
        with project.group(model_name):
            nodes[model_name] = GSCDB138Benchmark(model=model, model_name=model_name)
    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()


def test_gscdb138():
    build_project(repro=False)

