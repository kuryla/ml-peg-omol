"""Analyse GSCDB138 benchmark with per‑subset inspection."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import pytest

from ml_peg.analysis.utils.decorators import build_table
from ml_peg.analysis.utils.utils import mae
from ml_peg.app import APP_ROOT
from ml_peg.calcs import CALCS_ROOT
from ml_peg.models.get_models import get_model_names
from ml_peg.models.models import current_models


MODELS = get_model_names(current_models)
CALC_PATH = CALCS_ROOT / "molecular_properties" / "GSCDB138" / "outputs"
OUT_PATH = APP_ROOT / "data" / "molecular_properties" / "GSCDB138"


def _read_rows(model_dir: Path) -> list[dict]:
    """Read aggregated or per‑subset CSVs for one model."""
    rows: list[dict] = []
    agg = model_dir / "gscdb138_results.csv"
    if agg.exists():
        with open(agg) as fh:
            r = csv.DictReader(fh)
            for row in r:
                rows.append(row)
        return rows
    for p in sorted(model_dir.glob("*.csv")):
        with open(p) as fh:
            r = csv.DictReader(fh)
            for row in r:
                # allow subset implied by filename
                if not row.get("subset"):
                    row["subset"] = p.stem
                rows.append(row)
    return rows


def _collect() -> tuple[dict[str, dict[str, list[float]]], list[str]]:
    """
    Collect per‑subset reference and predictions across models.

    Returns
    -------
    (results, subsets)
      results[subset] = {"ref": [...], model: [...]}
    """
    results: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    subsets: set[str] = set()

    for model in MODELS:
        mdir = CALC_PATH / model
        if not mdir.exists():
            continue
        for row in _read_rows(mdir):
            try:
                subset = str(row["subset"]) if row.get("subset") else "unknown"
                subsets.add(subset)
                ref = float(row["ref"]) if row.get("ref") not in (None, "") else None
                pred = float(row["pred"]) if row.get("pred") not in (None, "") else None
                if ref is None or pred is None:
                    continue
                results[subset]["ref"].append(ref)
                results[subset][model].append(pred)
            except Exception:
                continue
    return results, sorted(subsets)


@pytest.fixture
def per_subset_results() -> dict[str, dict[str, list[float]]]:
    res, subsets = _collect()
    # Write manifest with per‑subset arrays for the app
    OUT_PATH.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH / "combined_results.json", "w") as fh:
        json.dump(res, fh)
    with open(OUT_PATH / "subsets.json", "w") as fh:
        json.dump(subsets, fh)

    # Also emit small metrics tables per subset for optional viewing
    tables_dir = OUT_PATH / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    for subset, data in res.items():
        metrics = {"MAE": {m: mae(data["ref"], data.get(m, [])) for m in MODELS}}
        # Build a minimal dash table JSON format
        rows = [{"MLIP": m, "MAE": metrics["MAE"][m], "id": m} for m in MODELS]
        table = {
            "data": rows,
            "columns": [{"name": i, "id": i} for i in ("MLIP", "MAE")],
            "tooltip_header": {"MLIP": "Model", "MAE": "Mean absolute error"},
        }
        with open(tables_dir / f"{subset}_metrics_table.json", "w") as fh:
            json.dump(table, fh)
    return res


@pytest.fixture
@build_table(
    filename=OUT_PATH / "gscdb138_metrics_table.json",
    metric_tooltips={
        "MLIP": "Model name",
        "MAE": "Mean absolute error across all subsets",
    },
)
def metrics(per_subset_results: dict[str, dict[str, list[float]]]) -> dict[str, dict]:
    # Overall MAE across concatenated subsets
    overall_ref: list[float] = []
    overall_pred: dict[str, list[float]] = {m: [] for m in MODELS}
    for subset, data in per_subset_results.items():
        overall_ref.extend(data["ref"])
        for m in MODELS:
            overall_pred[m].extend(data.get(m, []))
    overall = {m: mae(overall_ref, overall_pred[m]) if overall_pred[m] else None for m in MODELS}
    return {"MAE": overall}


def test_gscdb138(metrics: dict[str, dict]) -> None:
    return

