"""Run GSCDB138 app with per‑subset inspection."""

from __future__ import annotations

import json
from pathlib import Path

import plotly.graph_objects as go
from dash import Dash
from dash import Input, Output, callback, dcc
from dash.html import Div

from ml_peg.app import APP_ROOT
from ml_peg.app.base_app import BaseApp
from ml_peg.app.utils.load import rebuild_table
from ml_peg.models.get_models import get_model_names
from ml_peg.models.models import current_models


MODELS = get_model_names(current_models)
BENCHMARK_NAME = "GSCDB138"
DOCS_URL = "https://github.com/JiashuLiang/GSCDB138"
DATA_PATH = APP_ROOT / "data" / "molecular_properties" / "GSCDB138"


def _load_results() -> dict:
    p = DATA_PATH / "combined_results.json"
    if not p.exists():
        return {}
    with open(p) as fh:
        return json.load(fh)


def _load_subsets() -> list[str]:
    p = DATA_PATH / "subsets.json"
    if not p.exists():
        return []
    with open(p) as fh:
        return json.load(fh)


class GSCDB138App(BaseApp):
    """App for GSCDB138 with subset dropdown and per‑subset table/plot."""

    def register_callbacks(self) -> None:
        # Subset selector and containers already in layout
        @callback(
            Output(f"{BENCHMARK_NAME}-figure-placeholder", "children"),
            Output(f"{BENCHMARK_NAME}-subset-table-placeholder", "children"),
            Input(f"{BENCHMARK_NAME}-subset", "value"),
        )
        def update_subset(subset: str):
            res = _load_results()
            if not subset or subset not in res:
                return Div("Select a subset."), Div()

            data = res[subset]
            # build parity plot
            fig = go.Figure()
            ref = data.get("ref", [])
            for m in MODELS:
                pred = data.get(m, [])
                if pred:
                    fig.add_trace(
                        go.Scatter(x=pred, y=ref, mode="markers", name=m)
                    )
            if ref:
                x_min = min(min(data.get(m, [0])) for m in MODELS if data.get(m))
                x_max = max(max(data.get(m, [1])) for m in MODELS if data.get(m))
                r_min = min(ref)
                r_max = max(ref)
                lo = min(x_min, r_min)
                hi = max(x_max, r_max)
                fig.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines", showlegend=False))
            fig.update_layout(
                title=f"{subset} parity plot",
                xaxis_title="Predicted",
                yaxis_title="Reference",
            )

            # load subset table from prebuilt JSON
            table_path = DATA_PATH / "tables" / f"{subset}_metrics_table.json"
            if table_path.exists():
                table = rebuild_table(table_path, id=f"{BENCHMARK_NAME}-{subset}-table")
            else:
                table = Div("No metrics table available for this subset.")

            return dcc.Graph(figure=fig), table


def get_app() -> GSCDB138App:
    # Build overall metrics table
    metrics_table = rebuild_table(DATA_PATH / "gscdb138_metrics_table.json", id=f"{BENCHMARK_NAME}-table")
    # Compose layout extras
    subset_options = [
        {"label": s, "value": s} for s in _load_subsets()
    ]
    extras = [
        dcc.Dropdown(id=f"{BENCHMARK_NAME}-subset", options=subset_options, placeholder="Select a subset"),
        Div(id=f"{BENCHMARK_NAME}-figure-placeholder", className="graph-placeholder"),
        Div(id=f"{BENCHMARK_NAME}-subset-table-placeholder"),
    ]
    return GSCDB138App(
        name=BENCHMARK_NAME,
        description=(
            "Gold-Standard Chemical Database 138 with per‑subset inspection."
        ),
        docs_url=DOCS_URL,
        table_path=DATA_PATH / "gscdb138_metrics_table.json",
        extra_components=extras,
    )


if __name__ == "__main__":
    full_app = Dash(__name__, assets_folder=DATA_PATH.parent.parent)
    app = get_app()
    full_app.layout = app.layout
    app.register_callbacks()
    full_app.run(port=8071, debug=True)

