"""Analyse TMG145 RMSD results (bootstrap)."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ml_peg.analysis.utils.decorators import build_table
from ml_peg.app import APP_ROOT
from ml_peg.calcs import CALCS_ROOT
from ml_peg.models.get_models import get_model_names
from ml_peg.models.models import current_models


MODELS = get_model_names(current_models)
CALC_PATH = CALCS_ROOT / "tm_complexes" / "TMG145" / "outputs"
OUT_PATH = APP_ROOT / "data" / "tm_complexes" / "TMG145"


@pytest.fixture
@build_table(
    filename=OUT_PATH / "tmg145_metrics_table.json",
    metric_tooltips={
        "MLIP": "Model name",
        "Mean heavy RMSD": "Mean heavy-atom RMSD across TMG145 (Å)",
    },
)
def metrics() -> dict[str, dict]:
    results: dict[str, float] = {}
    for m in MODELS:
        p = CALC_PATH / m / "tmg145_results.csv"
        if not p.exists():
            results[m] = float("nan")
            continue
        vals=[]
        with open(p) as fh:
            r = csv.DictReader(fh)
            for row in r:
                try:
                    vals.append(float(row["heavy_rmsd"]))
                except Exception:
                    continue
        results[m] = sum(vals)/len(vals) if vals else float("nan")
    return {"Mean heavy RMSD": results}


def test_tmg145(metrics: dict[str, dict]) -> None:
    return

