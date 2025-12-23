"""Run L14 app."""

from __future__ import annotations

from pathlib import Path

from dash import Dash
from dash.html import Div

from ml_peg.app import APP_ROOT
from ml_peg.app.base_app import BaseApp
from ml_peg.app.utils.build_callbacks import (
    plot_from_table_column,
    struct_from_scatter,
)
from ml_peg.app.utils.load import read_plot
from ml_peg.models.get_models import get_model_names
from ml_peg.models.models import current_models

MODELS = get_model_names(current_models)
BENCHMARK_NAME = "L14"
DOCS_URL = (
    "https://ddmms.github.io/ml-peg/user_guide/benchmarks/supramolecular.html#L14"
)
DATA_PATH = APP_ROOT / "data" / "supramolecular" / "L14"


def _structure_paths(model_name: str) -> list[str]:
    """
    Build ordered list of structures for the scatter callback.

    Parameters
    ----------
    model_name
        Model name whose structures should be displayed.

    Returns
    -------
    list[str]
        Paths to XYZ files served through Dash assets.
    """
    struct_dir = DATA_PATH / model_name
    xyz_files = sorted(struct_dir.glob("*.xyz"))
    prefix = Path("assets") / "supramolecular" / "L14" / model_name
    return [str(prefix / path.name) for path in xyz_files]


class L14App(BaseApp):
    """L14 benchmark app layout and callbacks."""

    def register_callbacks(self) -> None:
        """Register callbacks to link metrics, parity plot, and structures."""
        scatter = read_plot(
            DATA_PATH / "figure_L14.json",
            id=f"{BENCHMARK_NAME}-figure",
        )

        # Structures are stored for each model under the benchmark data folder.
        structs = _structure_paths(MODELS[0])

        plot_from_table_column(
            table_id=self.table_id,
            plot_id=f"{BENCHMARK_NAME}-figure-placeholder",
            column_to_plot={
                "MAE": scatter,
                "RMSE": scatter,
            },
        )

        struct_from_scatter(
            scatter_id=f"{BENCHMARK_NAME}-figure",
            struct_id=f"{BENCHMARK_NAME}-struct-placeholder",
            structs=structs,
            mode="struct",
        )


def get_app() -> L14App:
    """
    Get L14 benchmark app layout and callback registration.

    Returns
    -------
    L14
        Benchmark layout and callback registration.
    """
    return L14App(
        name=BENCHMARK_NAME,
        description=(
            "Performance in predicting interaction energies for the "
            "L14 complexes that comprise the L14 benchmark."
        ),
        docs_url=DOCS_URL,
        table_path=DATA_PATH / "L14_metrics_table.json",
        extra_components=[
            Div(id=f"{BENCHMARK_NAME}-figure-placeholder"),
            Div(id=f"{BENCHMARK_NAME}-struct-placeholder"),
        ],
    )


if __name__ == "__main__":
    full_app = Dash(__name__, assets_folder=DATA_PATH.parent.parent)

    plf_app = get_app()
    full_app.layout = plf_app.layout
    plf_app.register_callbacks()

    full_app.run(port=8055, debug=True)
