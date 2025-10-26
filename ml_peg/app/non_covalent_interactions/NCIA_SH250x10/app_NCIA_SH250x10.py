"""Run NCIA_SH250x10 app (stub)."""

from __future__ import annotations

from dash import Dash
from dash.html import Div

from ml_peg.app import APP_ROOT
from ml_peg.app.base_app import BaseApp

BENCHMARK_NAME = "NCIA_SH250x10"
DOCS_URL = "https://ddmms.github.io/ml-peg/user_guide/benchmarks/non_covalent_interactions.html#ncia-sh250x10"
DATA_PATH = APP_ROOT / "data" / "non_covalent_interactions" / "NCIA_SH250x10"


class NCIA_SH250x10_App(BaseApp):
    def register_callbacks(self) -> None:
        return None


def get_app() -> NCIA_SH250x10_App:
    return NCIA_SH250x10_App(
        name=BENCHMARK_NAME,
        description=("NCIA SH250x10: sigma-hole interactions (table placeholder)."),
        docs_url=DOCS_URL,
        table_path=DATA_PATH / "ncia_sh250x10_metrics_table.json",
        extra_components=[Div()],
    )


if __name__ == "__main__":
    full_app = Dash(__name__, assets_folder=DATA_PATH.parent.parent)
    app = get_app()
    full_app.layout = app.layout
    app.register_callbacks()
    full_app.run(port=8069, debug=True)

