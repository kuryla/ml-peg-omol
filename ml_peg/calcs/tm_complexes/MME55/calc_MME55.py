"""Bootstrap calculations for the MME55 benchmark (metalloenzymes)."""

from __future__ import annotations

import csv
from pathlib import Path

import mlipx
import pytest
import zntrack

from ml_peg.calcs.utils.utils import chdir
from ml_peg.models.get_models import load_models
from ml_peg.models.models import current_models


MODELS = load_models(current_models)


def get_paths() -> tuple[Path, Path]:
    data = Path(__file__).parent / "data" / "results_MME55"
    out_base = Path(__file__).parent / "outputs"
    if not data.exists():
        raise FileNotFoundError("MME55 results not found under data/results_MME55")
    return data, out_base


class MME55Benchmark(zntrack.Node):
    model_name: str = zntrack.params()

    def run(self):
        data, out_base = get_paths()
        src = data / "all_detailed_results.csv"
        out_dir = out_base / self.model_name
        out_dir.mkdir(parents=True, exist_ok=True)
        dst = out_dir / "mme55_results.csv"

        # Normalize header to our expected columns
        with open(src) as fh, open(dst, "w", newline="") as fo:
            r = csv.DictReader(fh)
            w = csv.writer(fo)
            w.writerow(["reaction_id", "enzyme", "step_type", "ref_value", "calc_value", "error", "abs_error", "atoms_count"])
            for row in r:
                w.writerow([
                    row.get("reaction_id"),
                    row.get("enzyme"),
                    row.get("step_type"),
                    row.get("ref_value"),
                    row.get("calculated_value"),
                    row.get("error"),
                    row.get("abs_error"),
                    row.get("atoms_count"),
                ])


def build_project(repro: bool = False) -> None:
    project = mlipx.Project()
    nodes: dict[str, MME55Benchmark] = {}
    for model_name in MODELS:
        with project.group(model_name):
            nodes[model_name] = MME55Benchmark(model_name=model_name)
    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()
        for n in nodes.values():
            n.run()


def test_mme55():
    build_project(repro=False)

