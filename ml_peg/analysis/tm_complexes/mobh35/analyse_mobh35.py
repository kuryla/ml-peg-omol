"""Analyse MOBH35 reaction barrier benchmark."""

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
CALC_PATH = CALCS_ROOT / "single_point" / "mobh35" / "outputs"
OUT_PATH = APP_ROOT / "data" / "single_point" / "mobh35"

# Sort out the prefixes of reactants, products, ts
reactant_prefixes = {i: [f'r{i}'] for i in range(1, 36)}
product_prefixes = {i: [f'p{i}'] for i in range(1, 36)}
ts_prefixes = {i: [f'ts{i}'] for i in range(1, 36)}

for i in [1, 2, 32, 33]:
    reactant_prefixes[i] = [f'r{i}+.xyz']
    product_prefixes[i] = [f'p{i}+.xyz']
    ts_prefixes[i] = [f'ts{i}+.xyz']

reactant_prefixes[30] = ['r30_r31']
reactant_prefixes[31] = ['r30_r31']
reactant_prefixes[35] = ['p34_r35']
product_prefixes[34] = ['p34_r35']
reactant_prefixes[7] = ['p6_r7']
product_prefixes[6] = ['p6_r7']
reactant_prefixes[9] = ['p8_r9']
product_prefixes[8] = ['p8_r9']

product_prefixes[10].append('CO')
for i in [17, 18, 19, 20]:
    product_prefixes[i].append('PEt3')
product_prefixes[35].append('CH4')


def labels() -> list:
    """
    Get list of system names.

    Returns
    -------
    list
        List of all system names.
    """
    for model_name in MODELS:
        labels_list = [path.stem.replace('ts', '').replace('+', '') for path in sorted((CALC_PATH / model_name).glob("ts*"))]
        return labels_list


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_mobh35_barriers.json",
    title="Reaction barriers",
    x_label="Predicted barrier / eV",
    y_label="Reference barrier / eV",
    hoverdata={
        "Labels": labels(),
    },
)
def barrier_heights() -> dict[str, list]:
    """
    Get barrier heights for all systems.

    Returns
    -------
    dict[str, list]
        Dictionary of all reference and predicted barrier heights.
    """
    results = {"ref": []} | {mlip: [] for mlip in MODELS}
    ref_stored = False
    
    for model_name in MODELS:
        for label in labels():
            forward_barrier = 0
            reverse_barrier = 0

            # Write structures for app
            structs_dir = OUT_PATH / model_name
            structs_dir.mkdir(parents=True, exist_ok=True)

            for prefix in ts_prefixes:
                atoms = read(CALC_PATH / model_name / f'{prefix}.xyz')
                forward_barrier += atoms.info['model_energy']
                reverse_barrier += atoms.info['model_energy']
                if not ref_stored:
                    results['ref'].append(atoms.info['ref_forward_barrier'])
                    results['ref'].append(atoms.info['ref_reverse_barrier'])
                write(structs_dir / f"{prefix}.xyz", atoms)
            
            for prefix in reactant_prefixes:
                atoms = read(CALC_PATH / model_name / f'{prefix}.xyz')
                forward_barrier -= atoms.info['model_energy']
                write(structs_dir / f"{prefix}.xyz", atoms)
            
            for prefix in product_prefixes:
                atoms = read(CALC_PATH / model_name / f'{prefix}.xyz')
                reverse_barrier -= atoms.info['model_energy']
                write(structs_dir / f"{prefix}.xyz", atoms)

            results[model_name].append(forward_barrier)
            results[model_name].append(reverse_barrier)
        ref_stored = True
    return results


@pytest.fixture
def get_mae(barrier_heights) -> dict[str, float]:
    """
    Get mean absolute error for barrier heights.

    Parameters
    ----------
    barrier_heights
        Dictionary of reference and predicted barrier heights.

    Returns
    -------
    dict[str, float]
        Dictionary of predicted barrier height errors for all models.
    """
    results = {}
    for model_name in MODELS:
        results[model_name] = mae(
            barrier_heights["ref"], barrier_heights[model_name]
        )
    return results


@pytest.fixture
def get_rmse(barrier_heights) -> dict[str, float]:
    """
    Get root mean square error for barrier heights.

    Parameters
    ----------
    barrier_heights
        Dictionary of reference and predicted barrier heights.

    Returns
    -------
    dict[str, float]
        Dictionary of predicted barrier height errors for all models.
    """
    results = {}
    for model_name in MODELS:
        results[model_name] = rmse(
            barrier_heights["ref"], barrier_heights[model_name]
        )
    return results


@pytest.fixture
@build_table(
    filename=OUT_PATH / "mobh35_barriers_metrics_table.json",
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


def test_mobh35_barriers(metrics: dict[str, dict]) -> None:
    """
    Run bh9_barriers test.

    Parameters
    ----------
    metrics
        All new benchmark metric names and dictionary of values for each model.
    """
    return