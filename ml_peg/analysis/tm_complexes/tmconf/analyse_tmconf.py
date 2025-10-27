"""Analyse TMCONF benchmark. 10.1021/acs.jpca.2c02951"""

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


MODELS = load_models(current_models)

KCAL_TO_EV = units.kcal / units.mol
EV_TO_KCAL = 1 / KCAL_TO_EV
CALC_PATH = CALCS_ROOT / "single_point" / "tmconf" / "outputs"
OUT_PATH = APP_ROOT / "data" / "single_point" / "tmconf"
print(APP_ROOT)

def labels() -> list:
    """
    Get list of system names.

    Returns
    -------
    list
        List of all system names.
    """
    labels_list = [path.stem for path in sorted(CALC_PATH.glob('*.xyz'))]
    return labels_list


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_tmconf.json",
    title="Interaction energies",
    x_label="Predicted energy / eV",
    y_label="Reference energy / eV",
    hoverdata={
        "Labels": labels(),
    },
)
def interaction_energies() -> dict[str, list]:
    """
    Get interaction energies for all systems.

    Returns
    -------
    dict[str, list]
        Dictionary of all reference and predicted interaction energies.
    """
    results = {"ref": []} | {mlip: [] for mlip in MODELS}

    ref_stored = False

    for model_name in MODELS:
        for label in labels():
            atoms = read(CALC_PATH / model_name / f'{label}.xyz')
            # Get model and ref energies of the zero enegry conformation
            if label[-1] == '1':
                E_zero_model = atoms.info['model_energy']
                E_zero_ref = atoms.info['ref_energy']
                continue
            
            if not ref_stored:
                results['ref'].append(atoms.info['ref_energy'] - E_zero_ref)

            results[model_name].append(atoms.info['model_energy'] - E_zero_model)
            
            # Write structures for app
            structs_dir = OUT_PATH / model_name
            structs_dir.mkdir(parents=True, exist_ok=True)
            write(structs_dir / f"{label}.xyz", atoms)

        ref_stored = True
    return results


@pytest.fixture
def get_mae(energies) -> dict[str, float]:
    """
    Get mean absolute error for energies.

    Parameters
    ----------
    energies
        Dictionary of reference and predicted energies.

    Returns
    -------
    dict[str, float]
        Dictionary of predicted energy errors for all models.
    """
    results = {}
    for model_name in MODELS:
        results[model_name] = mae(
            energies["ref"], energies[model_name]
        )
    return results


@pytest.fixture
def get_rmse(energies) -> dict[str, float]:
    """
    Get root mean square error for energies.

    Parameters
    ----------
    energies
        Dictionary of reference and predicted energies.

    Returns
    -------
    dict[str, float]
        Dictionary of predicted energy errors for all models.
    """
    results = {}
    for model_name in MODELS:
        results[model_name] = rmse(
            energies["ref"], energies[model_name]
        )
    return results


@pytest.fixture
@build_table(
    filename=OUT_PATH / "tmconf_metrics_table.json",
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


def test_tmconf(metrics: dict[str, dict]) -> None:
    """
    Run tmconf test.

    Parameters
    ----------
    metrics
        All new benchmark metric names and dictionary of values for each model.
    """
    return