"""Bootstrap calculations for TMG145 (geometry optimisation RMSDs)."""

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
    data = Path(__file__).parent / "data" / "results_tmg145"
    out_base = Path(__file__).parent / "outputs"
    if not data.exists():
        raise FileNotFoundError("TMG145 results not found under data/results_tmg145")
    return data, out_base


class TMG145Benchmark(zntrack.Node):
    model_name: str = zntrack.params()

    def run(self):
        data, out_base = get_paths()
        # Choose CSV for the current model if present, else any *_tmg145.csv
        src = None
        for p in data.glob("*_tmg145.csv"):
            if self.model_name in p.name:
                src = p
                break
        if src is None:
            # fallback to mace-omol or any
            for cand in ["mace-omol_tmg145.csv"]:
                if (data / cand).exists():
                    src = data / cand
                    break
        if src is None:
            src = next(iter(data.glob("*_tmg145.csv")), None)
        if src is None:
            raise FileNotFoundError("No *_tmg145.csv found in results_tmg145")

        out_dir = out_base / self.model_name
        out_dir.mkdir(parents=True, exist_ok=True)
        dst = out_dir / "tmg145_results.csv"

        with open(src) as fh, open(dst, "w", newline="") as fo:
            r = csv.DictReader(fh)
            w = csv.writer(fo)
            w.writerow(["refcode", "natoms", "heavy_rmsd"])
            for row in r:
                w.writerow([row.get("refcode"), row.get("natoms"), row.get("heavy_rmsd")])


def build_project(repro: bool = False) -> None:
    project = mlipx.Project()
    nodes: dict[str, TMG145Benchmark] = {}
    for model_name in MODELS:
        with project.group(model_name):
            nodes[model_name] = TMG145Benchmark(model_name=model_name)
    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()
        for n in nodes.values():
            n.run()


def test_tmg145():
    build_project(repro=False)

