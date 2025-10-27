"""Run MOR41 app."""

from __future__ import annotations

from dash import Dash
from dash.html import Div

from ml_peg.app import APP_ROOT
from ml_peg.app.base_app import BaseApp
from ml_peg.app.utils.build_callbacks import (
    plot_from_table_column,
    struct_from_scatter,
)
from ml_peg.app.utils.load import read_plot
import json
from ml_peg.models.get_models import get_model_names
from ml_peg.models.models import current_models

MODELS = get_model_names(current_models)
BENCHMARK_NAME = "MOR41"
DOCS_URL = "https://ddmms.github.io/ml-peg/user_guide/benchmarks/tm_complexes.html#mor41"
DATA_PATH = APP_ROOT / "data" / "tm_complexes" / "MOR41"


class MOR41App(BaseApp):
    """MOR41 benchmark app layout and callbacks."""

    def register_callbacks(self) -> None:
        """Register callbacks to app."""
        scatter = read_plot(
            DATA_PATH / "figure_mor41.json",
            id=f"{BENCHMARK_NAME}-figure",
        )

        plot_from_table_column(
            table_id=self.table_id,
            plot_id=f"{BENCHMARK_NAME}-figure-placeholder",
            column_to_plot={"MAE": scatter},
        )

        # Load ids order to map clicks -> correct structure
        structs: list[str] = []
        try:
            with open(DATA_PATH / "ids.json") as fh:
                ids = json.load(fh)
            structs = [
                f"assets/tm_complexes/MOR41/{MODELS[0]}/{rid}.xyz" for rid in ids
            ]
        except Exception:
            pass

        struct_from_scatter(
            scatter_id=f"{BENCHMARK_NAME}-figure",
            struct_id=f"{BENCHMARK_NAME}-struct-placeholder",
            structs=structs,
            mode="struct",
        )


def get_app() -> MOR41App:
    """
    Get MOR41 benchmark app layout and callback registration.

    Returns
    -------
    MOR41App
        Benchmark layout and callback registration.
    """
    return MOR41App(
        name=BENCHMARK_NAME,
        description=(
            "Closed‑shell metal‑organic reaction energies (kcal/mol)."
        ),
        docs_url=DOCS_URL,
        table_path=DATA_PATH / "mor41_metrics_table.json",
        extra_components=[
            Div(id=f"{BENCHMARK_NAME}-figure-placeholder"),
            Div(id=f"{BENCHMARK_NAME}-struct-placeholder"),
        ],
    )


if __name__ == "__main__":
    full_app = Dash(__name__, assets_folder=DATA_PATH.parent.parent)
    app = get_app()
    full_app.layout = app.layout
    app.register_callbacks()
    full_app.run(port=8057, debug=True)
