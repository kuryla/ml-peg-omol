"""Analyse MME55 benchmark (bootstrap from provided results)."""

from __future__ import annotations

import csv
from pathlib import Path
from collections import defaultdict

import pytest

from ml_peg.analysis.utils.decorators import build_table, plot_parity
from ml_peg.analysis.utils.utils import mae
from ml_peg.app import APP_ROOT
from ml_peg.calcs import CALCS_ROOT
from ml_peg.models.get_models import get_model_names
from ml_peg.models.models import current_models


MODELS = get_model_names(current_models)
CALC_PATH = CALCS_ROOT / "tm_complexes" / "MME55" / "outputs"
OUT_PATH = APP_ROOT / "data" / "tm_complexes" / "MME55"


def reaction_ids() -> list[str]:
    for m in MODELS:
        p = CALC_PATH / m / "mme55_results.csv"
        if p.exists():
            with open(p) as fh:
                r = csv.DictReader(fh)
                return [row["reaction_id"] for row in r]
    return []


def _collect() -> dict[str, list]:
    results: dict[str, list] = {"ref": []} | {mlip: [] for mlip in MODELS}
    ref_stored = False
    for m in MODELS:
        p = CALC_PATH / m / "mme55_results.csv"
        if not p.exists():
            results[m] = []
            continue
        refs: list[float] = []
        preds: list[float] = []
        with open(p) as fh:
            r = csv.DictReader(fh)
            for row in r:
                refs.append(float(row["ref_value"]))
                preds.append(float(row["calc_value"]))
        results[m] = preds
        if not ref_stored:
            results["ref"] = refs
            ref_stored = True
    return results


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_mme55.json",
    title="MME55: reaction/barrier values",
    x_label="Predicted (kcal/mol)",
    y_label="Reference (kcal/mol)",
    hoverdata={"Reaction": reaction_ids()},
)
def mme55_values() -> dict[str, list]:
    return _collect()


@pytest.fixture
def mme55_mae(mme55_values: dict[str, list]) -> dict[str, float]:
    return {m: mae(mme55_values["ref"], mme55_values[m]) for m in MODELS}


@pytest.fixture
@build_table(
    filename=OUT_PATH / "mme55_metrics_table.json",
    metric_tooltips={
        "MLIP": "Model name",
        "MAE": "Mean absolute error across MME55 (kcal/mol)",
    },
)
def metrics(mme55_mae: dict[str, float]) -> dict[str, dict]:
    return {"MAE": mme55_mae}


def test_mme55(metrics: dict[str, dict]) -> None:
    return

