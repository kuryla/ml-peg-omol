"""Analyse QUID non‑covalent benchmark results."""

from __future__ import annotations

import csv
from pathlib import Path
import json

import pytest

from ml_peg.analysis.utils.decorators import build_table, plot_parity
from ml_peg.analysis.utils.utils import mae
from ml_peg.app import APP_ROOT
from ml_peg.calcs import CALCS_ROOT
from ml_peg.models.get_models import get_model_names
from ml_peg.models.models import current_models


MODELS = get_model_names(current_models)
CALC_PATH = CALCS_ROOT / "non_covalent_interactions" / "QUID" / "outputs"
OUT_PATH = APP_ROOT / "data" / "non_covalent_interactions" / "QUID"


def dimer_ids() -> list[str]:
    for m in MODELS:
        p = CALC_PATH / m / "quid_results.csv"
        if p.exists():
            with open(p) as fh:
                r = csv.DictReader(fh)
                return [row["id"] for row in r]
    return []


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_quid.json",
    title="QUID interaction energies",
    x_label="Predicted E_int (kcal/mol)",
    y_label="Reference E_int (kcal/mol)",
    hoverdata={"Dimer": dimer_ids()},
)
def energies() -> dict[str, list]:
    results: dict[str, list] = {"ref": []} | {mlip: [] for mlip in MODELS}
    ref_stored = False
    for m in MODELS:
        p = CALC_PATH / m / "quid_results.csv"
        if not p.exists():
            results[m] = []
            continue
        refs: list[float] = []
        preds: list[float] = []
        with open(p) as fh:
            r = csv.DictReader(fh)
            for row in r:
                refs.append(float(row["ref_eint_kcal"]))
                preds.append(float(row["calc_eint_kcal"]))
        results[m] = preds
        if not ref_stored:
            results["ref"] = refs
            ref_stored = True
    # Side-effect: materialize dimer structures for first model
    try:
        import h5py
        from ase import Atoms
        ids = dimer_ids()
        if ids:
            h5 = Path(__file__).parents[3] / "calcs" / "non_covalent_interactions" / "QUID" / "data" / "QUID.h5"
            dst_root = OUT_PATH / MODELS[0]
            dst_root.mkdir(parents=True, exist_ok=True)
            with h5py.File(h5, "r") as f:
                for name in ids:
                    g = f[name]
                    nums = g["atoms"]["dimer"][()]
                    pos = g["positions"]["dimer"][()]
                    Atoms(numbers=nums, positions=pos).write(dst_root / f"{name}.xyz")
            with open(OUT_PATH / "ids.json", "w") as fh:
                json.dump(ids, fh)
    except Exception:
        pass
    return results


@pytest.fixture
@build_table(
    filename=OUT_PATH / "quid_metrics_table.json",
    metric_tooltips={
        "MLIP": "Model name",
        "MAE": "Mean absolute error of E_int (kcal/mol)",
    },
)
def metrics(energies: dict[str, list]) -> dict[str, dict]:
    return {"MAE": {m: mae(energies["ref"], energies[m]) for m in MODELS}}


def test_quid(metrics: dict[str, dict]) -> None:
    return
