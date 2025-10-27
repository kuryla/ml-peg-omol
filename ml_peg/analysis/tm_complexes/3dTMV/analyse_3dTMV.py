"""Analyse 3dTMV benchmark results (ionization potentials)."""

from __future__ import annotations

import csv
from pathlib import Path
import json
from collections import defaultdict

import pytest

from ml_peg.analysis.utils.decorators import build_table, plot_parity
from ml_peg.analysis.utils.utils import mae
from ml_peg.app import APP_ROOT
from ml_peg.calcs import CALCS_ROOT
from ml_peg.models.get_models import get_model_names
from ml_peg.models.models import current_models


MODELS = get_model_names(current_models)
CALC_PATH = CALCS_ROOT / "tm_complexes" / "3dTMV" / "outputs"
OUT_PATH = APP_ROOT / "data" / "tm_complexes" / "3dTMV"


def complex_ids() -> list[int]:
    for model_name in MODELS:
        p = CALC_PATH / model_name / "tmv_results.csv"
        if p.exists():
            with open(p) as fh:
                r = csv.DictReader(fh)
                return [int(row["complex_id"]) for row in r]
    return []


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_ip.json",
    title="Vertical ionization potentials (kcal/mol)",
    x_label="Predicted IP (kcal/mol)",
    y_label="Reference IP (kcal/mol)",
    hoverdata={
        "Complex": complex_ids(),
    },
)
def ips() -> dict[str, list]:
    results: dict[str, list] = {"ref": []} | {mlip: [] for mlip in MODELS}
    ref_stored = False
    for model_name in MODELS:
        p = CALC_PATH / model_name / "tmv_results.csv"
        if not p.exists():
            results[model_name] = []
            continue
        refs: list[float] = []
        pred: list[float] = []
        with open(p) as fh:
            r = csv.DictReader(fh)
            for row in r:
                refs.append(float(row["IP_ref"]))
                pred.append(float(row["IP_model"]))
        results[model_name] = pred
        if not ref_stored:
            results["ref"] = refs
            ref_stored = True

    # Side-effect: copy structures for first available model into app assets dir
    try:
        # Determine ID order used in plot/table
        ids = complex_ids()
        if ids:
            src_root = Path(__file__).parents[3] / "calcs" / "tm_complexes" / "3dTMV" / "data" / "3dtmv_structures"
            dst_root = OUT_PATH / MODELS[0]
            dst_root.mkdir(parents=True, exist_ok=True)
            for cid in ids:
                src = src_root / str(cid) / "struc.xyz"
                if src.exists():
                    # Write a copy under assets data
                    from ase.io import read, write

                    atoms = read(src)
                    write(dst_root / f"{cid}.xyz", atoms)
            # Save the order to preserve index mapping in the app
            with open(OUT_PATH / "ids.json", "w") as fh:
                json.dump(ids, fh)
    except Exception:
        pass
    return results


def _read_by_subset(model_name: str) -> dict[str, list[float]]:
    p = CALC_PATH / model_name / "tmv_results.csv"
    by_subset: dict[str, list[tuple[float, float]]] = defaultdict(list)
    if not p.exists():
        return {}
    with open(p) as fh:
        r = csv.DictReader(fh)
        for row in r:
            by_subset[row["subset"]].append((float(row["IP_ref"]), float(row["IP_model"])))
    # Convert to mae per subset later
    return by_subset


@pytest.fixture
def tmv_mae(ips: dict[str, list]) -> dict[str, float]:
    return {m: mae(ips["ref"], ips[m]) for m in MODELS}


@pytest.fixture
def tmv_mae_sr() -> dict[str, float]:
    res = {}
    for model_name in MODELS:
        by_subset = _read_by_subset(model_name)
        if "SR" in by_subset:
            ref, pred = zip(*by_subset["SR"], strict=True)
            res[model_name] = mae(list(ref), list(pred))
        else:
            res[model_name] = float("nan")
    return res


@pytest.fixture
def tmv_mae_srmr() -> dict[str, float]:
    res = {}
    for model_name in MODELS:
        by_subset = _read_by_subset(model_name)
        if "SR/MR" in by_subset:
            ref, pred = zip(*by_subset["SR/MR"], strict=True)
            res[model_name] = mae(list(ref), list(pred))
        else:
            res[model_name] = float("nan")
    return res


@pytest.fixture
def tmv_mae_mr() -> dict[str, float]:
    res = {}
    for model_name in MODELS:
        by_subset = _read_by_subset(model_name)
        if "MR" in by_subset:
            ref, pred = zip(*by_subset["MR"], strict=True)
            res[model_name] = mae(list(ref), list(pred))
        else:
            res[model_name] = float("nan")
    return res


@pytest.fixture
@build_table(
    filename=OUT_PATH / "3dtmv_metrics_table.json",
    metric_tooltips={
        "MLIP": "Model name",
        "MAE (All)": "Mean absolute error across all complexes (kcal/mol)",
        "MAE (SR)": "MAE over SR subset",
        "MAE (SR/MR)": "MAE over SR/MR subset",
        "MAE (MR)": "MAE over MR subset",
    },
)
def metrics(
    tmv_mae: dict[str, float], tmv_mae_sr: dict[str, float], tmv_mae_srmr: dict[str, float], tmv_mae_mr: dict[str, float]
) -> dict[str, dict]:
    return {
        "MAE (All)": tmv_mae,
        "MAE (SR)": tmv_mae_sr,
        "MAE (SR/MR)": tmv_mae_srmr,
        "MAE (MR)": tmv_mae_mr,
    }


def test_3dtmv(metrics: dict[str, dict]) -> None:
    return
