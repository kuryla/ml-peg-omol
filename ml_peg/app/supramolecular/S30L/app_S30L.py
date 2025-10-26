"""Run S30L app (stub)."""

from __future__ import annotations

from dash import Dash
from dash.html import Div

from ml_peg.app import APP_ROOT
from ml_peg.app.base_app import BaseApp
from ml_peg.models.get_models import get_model_names
from ml_peg.models.models import current_models

MODELS = get_model_names(current_models)
BENCHMARK_NAME = "S30L"
DOCS_URL = "https://ddmms.github.io/ml-peg/user_guide/benchmarks/supramolecular.html#s30l"
DATA_PATH = APP_ROOT / "data" / "supramolecular" / "S30L"


class S30LApp(BaseApp):
    """S30L benchmark app layout and callbacks (table-only for now)."""

    def register_callbacks(self) -> None:  # no-op until analysis is added
        return None


def get_app() -> S30LApp:
    """
    Get S30L benchmark app layout.
    """
    return S30LApp(
        name=BENCHMARK_NAME,
        description=(
            "Supramolecular complexes interaction energies (table placeholder)."
        ),
        docs_url=DOCS_URL,
        table_path=DATA_PATH / "s30l_metrics_table.json",
        extra_components=[Div()],
    )


if __name__ == "__main__":
    full_app = Dash(__name__, assets_folder=DATA_PATH.parent.parent)
    app = get_app()
    full_app.layout = app.layout
    app.register_callbacks()
    full_app.run(port=8062, debug=True)

