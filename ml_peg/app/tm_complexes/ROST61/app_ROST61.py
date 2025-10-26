"""Run ROST61 app."""

from __future__ import annotations

from dash import Dash
from dash.html import Div

from ml_peg.app import APP_ROOT
from ml_peg.app.base_app import BaseApp
from ml_peg.app.utils.build_callbacks import (
    plot_from_table_column,
)
from ml_peg.app.utils.load import read_plot
from ml_peg.models.get_models import get_model_names
from ml_peg.models.models import current_models

MODELS = get_model_names(current_models)
BENCHMARK_NAME = "ROST61"
DOCS_URL = "https://ddmms.github.io/ml-peg/user_guide/benchmarks/tm_complexes.html#rost61"
DATA_PATH = APP_ROOT / "data" / "tm_complexes" / "ROST61"


class ROST61App(BaseApp):
    """ROST61 benchmark app layout and callbacks."""

    def register_callbacks(self) -> None:
        """Register callbacks to app."""
        scatter = read_plot(
            DATA_PATH / "figure_reaction_energies.json",
            id=f"{BENCHMARK_NAME}-figure",
        )

        plot_from_table_column(
            table_id=self.table_id,
            plot_id=f"{BENCHMARK_NAME}-figure-placeholder",
            column_to_plot={"MAE": scatter},
        )


def get_app() -> ROST61App:
    """
    Get ROST61 benchmark app layout and callback registration.

    Returns
    -------
    ROST61App
        Benchmark layout and callback registration.
    """
    return ROST61App(
        name=BENCHMARK_NAME,
        description=(
            "Open‑shell organometallic reaction energies (kcal/mol)."
        ),
        docs_url=DOCS_URL,
        table_path=DATA_PATH / "rost61_metrics_table.json",
        extra_components=[
            Div(id=f"{BENCHMARK_NAME}-figure-placeholder"),
        ],
    )


if __name__ == "__main__":
    full_app = Dash(__name__, assets_folder=DATA_PATH.parent.parent)
    app = get_app()
    full_app.layout = app.layout
    app.register_callbacks()
    full_app.run(port=8058, debug=True)

