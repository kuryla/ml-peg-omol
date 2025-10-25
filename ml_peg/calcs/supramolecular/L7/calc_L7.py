"""Bootstrap calculations for L7 (supramolecular interactions)."""

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
    data = Path(__file__).parent / "data" / "results_l7"
    out_base = Path(__file__).parent / "outputs"
    if not data.exists():
        raise FileNotFoundError("L7 results not found under data/results_l7")
    return data, out_base


class L7Benchmark(zntrack.Node):
    model_name: str = zntrack.params()

    def run(self):
        data, out_base = get_paths()
        # Prefer a file that contains ref and model values in kcal
        src = None
        for p in data.glob("*.csv"):
            # pick a file with columns E_int_ref_kcal and E_int_model_kcal if present
            with open(p) as fh:
                header = fh.readline()
                if "E_int_ref_kcal" in header and "E_int_model_kcal" in header:
                    src = p
                    break
        if src is None:
            raise FileNotFoundError("Could not find a results CSV with ref/model kcal columns in results_l7")

        out_dir = out_base / self.model_name
        out_dir.mkdir(parents=True, exist_ok=True)
        dst = out_dir / "l7_results.csv"

        with open(src) as fh, open(dst, "w", newline="") as fo:
            r = csv.DictReader(fh)
            w = csv.writer(fo)
            w.writerow(["system", "ref_kcal", "calc_kcal", "error_kcal", "abs_error_kcal"])
            for row in r:
                if row.get("system") == "Model":  # skip summary rows
                    continue
                refk = float(row["E_int_ref_kcal"]) if row.get("E_int_ref_kcal") else float(row["E_int_ref_eV"]) * 23.0605419453
                calck = float(row["E_int_model_kcal"]) if row.get("E_int_model_kcal") else float(row["E_int_model_eV"]) * 23.0605419453
                errk = calck - refk
                w.writerow([row.get("system"), refk, calck, errk, abs(errk)])


def build_project(repro: bool = False) -> None:
    project = mlipx.Project()
    nodes: dict[str, L7Benchmark] = {}
    for model_name in MODELS:
        with project.group(model_name):
            nodes[model_name] = L7Benchmark(model_name=model_name)
    if repro:
        with chdir(Path(__file__).parent):
            project.repro(build=True, force=True)
    else:
        project.build()
        for n in nodes.values():
            n.run()


def test_l7():
    build_project(repro=False)

