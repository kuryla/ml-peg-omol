"""Analyse ROST61 benchmark outputs (tm_complexes category)."""

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
CALC_PATH = CALCS_ROOT / "tm_complexes" / "ROST61" / "outputs"
OUT_PATH = APP_ROOT / "data" / "tm_complexes" / "ROST61"


def reaction_ids() -> list[str]:
    for model_name in MODELS:
        path = CALC_PATH / model_name / "rost61_results.csv"
        if path.exists():
            with open(path) as fh:
                reader = csv.DictReader(fh)
                return [row["reaction_id"] for row in reader]
    return []


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_reaction_energies.json",
    title="ROST61 Reaction Energies",
    x_label="Predicted ΔE (kcal/mol)",
    y_label="Reference ΔE (kcal/mol)",
    hoverdata={
        "Reaction": reaction_ids(),
    },
)
def reaction_energies() -> dict[str, list]:
    results: dict[str, list] = {"ref": []} | {mlip: [] for mlip in MODELS}
    ref_stored = False
    for model_name in MODELS:
        path = CALC_PATH / model_name / "rost61_results.csv"
        if not path.exists():
            results[model_name] = []
            continue
        ref_vals: list[float] = []
        model_vals: list[float] = []
        with open(path) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                ref_vals.append(float(row["ref_energy"]))
                model_vals.append(float(row["calc_energy"]))
        results[model_name] = model_vals
        if not ref_stored:
            results["ref"] = ref_vals
            ref_stored = True
    return results


@pytest.fixture
def rost61_mae(reaction_energies: dict[str, list]) -> dict[str, float]:
    return {m: mae(reaction_energies["ref"], reaction_energies[m]) for m in MODELS}


@pytest.fixture
@build_table(
    filename=OUT_PATH / "rost61_metrics_table.json",
    metric_tooltips={
        "MLIP": "Name of the model",
        "MAE": "Mean Absolute Error (kcal/mol)",
    },
)
def metrics(rost61_mae: dict[str, float]) -> dict[str, dict]:
    return {"MAE": rost61_mae}


def test_rost61(metrics: dict[str, dict]) -> None:
    return

