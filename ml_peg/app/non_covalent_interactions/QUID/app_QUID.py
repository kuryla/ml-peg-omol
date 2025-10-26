"""Run QUID app."""

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
BENCHMARK_NAME = "QUID"
DOCS_URL = "https://ddmms.github.io/ml-peg/user_guide/benchmarks/non_covalent_interactions.html#quid"
DATA_PATH = APP_ROOT / "data" / "non_covalent_interactions" / "QUID"


class QUIDApp(BaseApp):
    """QUID benchmark app layout and callbacks."""

    def register_callbacks(self) -> None:
        """Register callbacks to app."""
        scatter = read_plot(
            DATA_PATH / "figure_quid.json",
            id=f"{BENCHMARK_NAME}-figure",
        )

        plot_from_table_column(
            table_id=self.table_id,
            plot_id=f"{BENCHMARK_NAME}-figure-placeholder",
            column_to_plot={"MAE": scatter},
        )


def get_app() -> QUIDApp:
    """
    Get QUID benchmark app layout and callback registration.

    Returns
    -------
    QUIDApp
        Benchmark layout and callback registration.
    """
    return QUIDApp(
        name=BENCHMARK_NAME,
        description=(
            "Interaction energies for a diverse set of dimers (kcal/mol)."
        ),
        docs_url=DOCS_URL,
        table_path=DATA_PATH / "quid_metrics_table.json",
        extra_components=[
            Div(id=f"{BENCHMARK_NAME}-figure-placeholder"),
        ],
    )


if __name__ == "__main__":
    full_app = Dash(__name__, assets_folder=DATA_PATH.parent.parent)
    app = get_app()
    full_app.layout = app.layout
    app.register_callbacks()
    full_app.run(port=8059, debug=True)

