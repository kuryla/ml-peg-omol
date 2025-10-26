"""Run NCIA_R739x5 app (stub)."""

from __future__ import annotations

from dash import Dash
from dash.html import Div

from ml_peg.app import APP_ROOT
from ml_peg.app.base_app import BaseApp

BENCHMARK_NAME = "NCIA_R739x5"
DOCS_URL = "https://ddmms.github.io/ml-peg/user_guide/benchmarks/non_covalent_interactions.html#ncia-r739x5"
DATA_PATH = APP_ROOT / "data" / "non_covalent_interactions" / "NCIA_R739x5"


class NCIA_R739x5_App(BaseApp):
    def register_callbacks(self) -> None:
        return None


def get_app() -> NCIA_R739x5_App:
    return NCIA_R739x5_App(
        name=BENCHMARK_NAME,
        description=("NCIA R739x5: repulsive contacts (table placeholder)."),
        docs_url=DOCS_URL,
        table_path=DATA_PATH / "ncia_r739x5_metrics_table.json",
        extra_components=[Div()],
    )


if __name__ == "__main__":
    full_app = Dash(__name__, assets_folder=DATA_PATH.parent.parent)
    app = get_app()
    full_app.layout = app.layout
    app.register_callbacks()
    full_app.run(port=8068, debug=True)

