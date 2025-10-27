"""Analyse MOR41 benchmark results."""

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
CALC_PATH = CALCS_ROOT / "tm_complexes" / "MOR41" / "outputs"
OUT_PATH = APP_ROOT / "data" / "tm_complexes" / "MOR41"


def reaction_numbers() -> list[int]:
    for model_name in MODELS:
        p = CALC_PATH / model_name / "mor41_results.csv"
        if p.exists():
            with open(p) as fh:
                r = csv.DictReader(fh)
                return [int(row["reaction_number"]) for row in r]
    return []


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_mor41.json",
    title="MOR41 reaction energies",
    x_label="Predicted ΔE (kcal/mol)",
    y_label="Reference ΔE (kcal/mol)",
    hoverdata={"Reaction": reaction_numbers()},
)
def energies() -> dict[str, list]:
    results: dict[str, list] = {"ref": []} | {mlip: [] for mlip in MODELS}
    ref_stored = False
    for model_name in MODELS:
        p = CALC_PATH / model_name / "mor41_results.csv"
        if not p.exists():
            results[model_name] = []
            continue
        refs: list[float] = []
        preds: list[float] = []
        with open(p) as fh:
            r = csv.DictReader(fh)
            for row in r:
                refs.append(float(row["ref_energy"]))
                preds.append(float(row["calc_energy"]))
        results[model_name] = preds
        if not ref_stored:
            results["ref"] = refs
            ref_stored = True
    # Side-effect: copy representative product structure for each reaction
    try:
        from ase.io import read, write
        from ml_peg.calcs.tm_complexes.MOR41.calc_MOR41 import create_reaction_mapping

        nums = reaction_numbers()
        if nums:
            mapping = create_reaction_mapping()
            src_root = Path(__file__).parents[3] / "calcs" / "tm_complexes" / "MOR41" / "data" / "geometries"
            dst_root = OUT_PATH / MODELS[0]
            dst_root.mkdir(parents=True, exist_ok=True)
            ids_order = []
            for n in nums:
                rid = f"MOR41_{n:02d}"
                ids_order.append(rid)
                info = mapping.get(rid, {})
                # Prefer first product molecule as representative
                products = info.get("products", [])
                mol = products[0][0] if products else None
                if mol:
                    f = src_root / mol / "mol.xyz"
                    if f.exists():
                        at = read(f)
                        write(dst_root / f"{rid}.xyz", at)
            with open(OUT_PATH / "ids.json", "w") as fh:
                json.dump(ids_order, fh)
    except Exception:
        pass
    return results


def _by_type(model_name: str) -> dict[str, list[tuple[float, float]]]:
    p = CALC_PATH / model_name / "mor41_results.csv"
    d: dict[str, list[tuple[float, float]]] = defaultdict(list)
    if not p.exists():
        return {}
    with open(p) as fh:
        r = csv.DictReader(fh)
        for row in r:
            d[row["reaction_type"]].append((float(row["ref_energy"]), float(row["calc_energy"])) )
    return d


@pytest.fixture
def mor41_mae(energies: dict[str, list]) -> dict[str, float]:
    return {m: mae(energies["ref"], energies[m]) for m in MODELS}


@pytest.fixture
def mor41_mae_oxidative() -> dict[str, float]:
    res={}
    for m in MODELS:
        d=_by_type(m)
        if "oxidative addition" in d:
            ref,pred=zip(*d["oxidative addition"], strict=True)
            res[m]=mae(list(ref), list(pred))
        else:
            res[m]=float("nan")
    return res


@pytest.fixture
@build_table(
    filename=OUT_PATH / "mor41_metrics_table.json",
    metric_tooltips={
        "MLIP": "Model name",
        "MAE": "Mean absolute error across MOR41 (kcal/mol)",
        "MAE (Oxidative addition)": "MAE for oxidative addition subset",
    },
)
def metrics(mor41_mae: dict[str, float], mor41_mae_oxidative: dict[str, float]) -> dict[str, dict]:
    return {"MAE": mor41_mae, "MAE (Oxidative addition)": mor41_mae_oxidative}


def test_mor41(metrics: dict[str, dict]) -> None:
    return
