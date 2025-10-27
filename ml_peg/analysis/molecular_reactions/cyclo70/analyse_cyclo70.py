"""Analyse CYCLO70 reaction barrier benchmark.

10.1021/acs.jctc.5c00925

"""

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

MODELS = load_models(current_models)

KCAL_TO_EV = units.kcal / units.mol
EV_TO_KCAL = 1 / KCAL_TO_EV
CALC_PATH = CALCS_ROOT / "single_point" / "cyclo70" / "outputs"
OUT_PATH = APP_ROOT / "data" / "single_point" / "cyclo70"
print(APP_ROOT)

def labels() -> list:
    """
    Get list of system names.

    Returns
    -------
    list
        List of all system names.
    """
    for model_name in MODELS:
        labels_list = [path.stem for path in sorted((CALC_PATH / model_name).glob("TS_*.xyz"))]
        return labels_list


@pytest.fixture
@plot_parity(
    filename=OUT_PATH / "figure_cyclo70_barriers.json",
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
            atoms_ts = read(CALC_PATH / model_name / f'{label}.xyz')

            results[model_name].append(atoms_ts.info['model_forward_bh'])
            results[model_name].append(atoms_ts.info['model_reverse_bh'])
            
            if not ref_stored:
                results['ref'].append(atoms_ts.info['ref_forward_bh'])
                results['ref'].append(atoms_ts.info['ref_reverse_bh'])
            
            # Write structures for app
            structs_dir = OUT_PATH / model_name
            structs_dir.mkdir(parents=True, exist_ok=True)
            write(structs_dir / f"{label}.xyz", atoms_ts)
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
    filename=OUT_PATH / "cyclo70_barriers_metrics_table.json",
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


def test_cyclo70_barriers(metrics: dict[str, dict]) -> None:
    """
    Run cyclo70 barriers test.

    Parameters
    ----------
    metrics
        All new benchmark metric names and dictionary of values for each model.
    """
    return