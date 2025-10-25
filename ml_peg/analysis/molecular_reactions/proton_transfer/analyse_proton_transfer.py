"""Analyse biochemical proton‑transfer results."""

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
CALC_PATH = CALCS_ROOT / "molecular_reactions" / "proton_transfer" / "outputs"
OUT_PATH = APP_ROOT / "data" / "molecular_reactions" / "proton_transfer"


def case_ids() -> list[str]:
    for m in MODELS:
        p = CALC_PATH / m / "pt_results.csv"
        if p.exists():
            with open(p) as fh:
                r = csv.DictReader(fh)
                return [row["case_id"] for row in r]
    return []


def _collect(metric: str) -> dict[str, list]:
    results: dict[str, list] = {"ref": []} | {mlip: [] for mlip in MODELS}
    ref_stored = False
    for m in MODELS:
        p = CALC_PATH / m / "pt_results.csv"
        if not p.exists():
            results[m] = []
            continue
        refs: list[float] = []
        preds: list[float] = []
        with open(p) as fh:
            r = csv.DictReader(fh)
            for row in r:
                ref = row["ref_rxn"] if metric == "rxn" else row["ref_barrier"]
                pred = row["calc_rxn"] if metric == "rxn" else row["calc_barrier"]
                if ref == "" or ref == "None":
                    # skip if no reference
                    continue
                refs.append(float(ref))
                preds.append(float(pred))
        results[m] = preds
        if not ref_stored:
            results["ref"] = refs
            ref_stored = True
    return results


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_pt_rxn.json",
    title="Proton‑transfer reaction energies",
    x_label="Predicted ΔE (kcal/mol)",
    y_label="Reference ΔE (kcal/mol)",
    hoverdata={"Case": case_ids()},
)
def pt_rxn() -> dict[str, list]:
    return _collect("rxn")


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_pt_barrier.json",
    title="Proton‑transfer forward barrier",
    x_label="Predicted barrier (kcal/mol)",
    y_label="Reference barrier (kcal/mol)",
    hoverdata={"Case": case_ids()},
)
def pt_barrier() -> dict[str, list]:
    return _collect("barrier")


@pytest.fixture
def pt_mae_rxn(pt_rxn: dict[str, list]) -> dict[str, float]:
    return {m: mae(pt_rxn["ref"], pt_rxn[m]) for m in MODELS}


@pytest.fixture
def pt_mae_barrier(pt_barrier: dict[str, list]) -> dict[str, float]:
    return {m: mae(pt_barrier["ref"], pt_barrier[m]) for m in MODELS}


@pytest.fixture
@build_table(
    filename=OUT_PATH / "pt_metrics_table.json",
    metric_tooltips={
        "MLIP": "Model name",
        "MAE (rxn)": "Mean absolute error for reaction energies (kcal/mol)",
        "MAE (barrier)": "Mean absolute error for barriers (kcal/mol)",
    },
)
def metrics(pt_mae_rxn: dict[str, float], pt_mae_barrier: dict[str, float]) -> dict[str, dict]:
    return {"MAE (rxn)": pt_mae_rxn, "MAE (barrier)": pt_mae_barrier}


def test_proton_transfer(metrics: dict[str, dict]) -> None:
    return

