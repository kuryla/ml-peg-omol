"""Analyse L7 results."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ml_peg.analysis.utils.decorators import build_table, plot_parity
from ml_peg.analysis.utils.utils import mae
from ml_peg.app import APP_ROOT
from ml_peg.calcs import CALCS_ROOT
from ml_peg.models.get_models import get_model_names
from ml_peg.models.models import current_models


MODELS = get_model_names(current_models)
CALC_PATH = CALCS_ROOT / "supramolecular" / "L7" / "outputs"
OUT_PATH = APP_ROOT / "data" / "supramolecular" / "L7"


def system_ids() -> list[str]:
    for m in MODELS:
        p = CALC_PATH / m / "l7_results.csv"
        if p.exists():
            with open(p) as fh:
                r = csv.DictReader(fh)
                return [row["system"] for row in r]
    return []


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_l7.json",
    title="L7 interaction energies",
    x_label="Predicted E_int (kcal/mol)",
    y_label="Reference E_int (kcal/mol)",
    hoverdata={"System": system_ids()},
)
def energies() -> dict[str, list]:
    results: dict[str, list] = {"ref": []} | {mlip: [] for mlip in MODELS}
    ref_stored = False
    for m in MODELS:
        p = CALC_PATH / m / "l7_results.csv"
        if not p.exists():
            results[m] = []
            continue
        refs: list[float] = []
        preds: list[float] = []
        with open(p) as fh:
            r = csv.DictReader(fh)
            for row in r:
                refs.append(float(row["ref_kcal"]))
                preds.append(float(row["calc_kcal"]))
        results[m] = preds
        if not ref_stored:
            results["ref"] = refs
            ref_stored = True
    return results


@pytest.fixture
@build_table(
    filename=OUT_PATH / "l7_metrics_table.json",
    metric_tooltips={"MLIP": "Model name", "MAE": "Mean absolute error (kcal/mol)"},
)
def metrics(energies: dict[str, list]) -> dict[str, dict]:
    return {"MAE": {m: mae(energies["ref"], energies[m]) for m in MODELS}}


def test_l7(metrics: dict[str, dict]) -> None:
    return

