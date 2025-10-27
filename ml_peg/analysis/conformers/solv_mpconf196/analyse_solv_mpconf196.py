"""Analyse solv-MPCONF196 conformer benchmark. 10.1002/jcc.27248"""

from __future__ import annotations

from ase.io import read, write
import pytest

from ml_peg.analysis.utils.decorators import build_table, plot_parity
from ml_peg.analysis.utils.utils import mae, rmse
from ml_peg.app import APP_ROOT
from ml_peg.calcs import CALCS_ROOT
from ml_peg.models.get_models import load_models
from ml_peg.models.models import current_models
from ml_peg.calcs.utils.utils import chdir, get_benchmark_data
from ase import units
import json
from collections import defaultdict
import numpy as np

MODELS = load_models(current_models)

KCAL_TO_EV = units.kcal / units.mol
EV_TO_KCAL = 1 / KCAL_TO_EV
CALC_PATH = CALCS_ROOT / "single_point" / "solv_mpconf196" / "outputs"
OUT_PATH = APP_ROOT / "data" / "single_point" / "solv_mpconf196"
MOLECULES = ['FGG', 'GFA', 'GGF', 'WG', 'WGG', 'CAMVES', 'CHPSAR', 'COHVAW', 'GS464992', 'GS557577', 'POXTRD', 'SANGLI', 'YIVNOG']


def labels() -> list:
    """
    Get list of system names.

    Returns
    -------
    list
        List of all system names.
    """
    for model_name in MODELS:
        labels_list = [path.stem for path in sorted((CALC_PATH / model_name).glob("*.xyz"))]
        return labels_list


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_solv_mpconf196.json",
    title="Energies",
    x_label="Predicted energy / eV",
    y_label="Reference energy / eV",
    hoverdata={
        "Labels": labels(),
    },
)
def conformer_energies() -> dict[str, list]:
    """
    Get conformer energies for all systems.

    Returns
    -------
    dict[str, list]
        Dictionary of all reference and predicted energies.
    """
    results = {"ref": []} | {mlip: [] for mlip in MODELS}
    ref_stored = False

    for model_name in MODELS:

        for label in labels():
            atoms = read(CALC_PATH / model_name / f'{label}.xyz')
            results[model_name].append(atoms.info['model_rel_energy'])
            if not ref_stored:
                results['ref'].append(atoms.info['ref_energy'])
            
            # Write structures for app
            structs_dir = OUT_PATH / model_name
            structs_dir.mkdir(parents=True, exist_ok=True)
            write(structs_dir / f"{label}.xyz", atoms)
        ref_stored = True
    return results


@pytest.fixture
def get_mae(conformer_energies) -> dict[str, float]:
    """
    Get mean absolute error for conformer energies.

    Parameters
    ----------
    conformer_energies
        Dictionary of reference and predicted conformer energies.

    Returns
    -------
    dict[str, float]
        Dictionary of predicted conformer energies errors for all models.
    """
    results = {}
    for model_name in MODELS:
        results[model_name] = mae(
            conformer_energies["ref"], conformer_energies[model_name]
        )
    return results


@pytest.fixture
def get_rmse(conformer_energies) -> dict[str, float]:
    """
    Get root mean square error for conformer energies.

    Parameters
    ----------
    conformer energies
        Dictionary of reference and predicted conformer energies.

    Returns
    -------
    dict[str, float]
        Dictionary of predicted conformer energies errors for all models.
    """
    results = {}
    for model_name in MODELS:
        results[model_name] = rmse(
            conformer_energies["ref"], conformer_energies[model_name]
        )
    return results


@pytest.fixture
@build_table(
    filename=OUT_PATH / "solv_mpconf196_metrics_table.json",
    metric_tooltips={
        "Model": "Name of the model",
        "MAE": "Mean Absolute Error (eV)",
        "RMSE": "Root Mean Square Error (eV)",
    },
)
def metrics(get_mae: dict[str, float], get_rmse: dict[str, float]) -> dict[str, dict]:
    """
    Get all metrics.

    Parameters
    ----------
    mae
        Mean absolute errors for all models.
    
    rmse
        Root Mean Square Error for all models

    Returns
    -------
    dict[str, dict]
        Metric names and values for all models.
    """
    return {
        "MAE/eV": get_mae,
        "RMSE/eV": get_rmse,
    }


def test_solv_mpconf196(metrics: dict[str, dict]) -> None:
    """
    Run solvMPCONF196 test.

    Parameters
    ----------
    metrics
        All new benchmark metric names and dictionary of values for each model.
    """
    return