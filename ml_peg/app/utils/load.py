"""Helpers to load components into Dash app."""

from __future__ import annotations

import json
from pathlib import Path

from dash.dash_table import DataTable
from dash.dcc import Graph
from plotly.io import read_json

<<<<<<< Updated upstream
from ml_peg.analysis.utils.utils import calc_metric_scores, get_table_style
from ml_peg.app.utils.utils import (
    calculate_column_widths,
    clean_thresholds,
    is_numeric_column,
    rank_format,
    sig_fig_format,
)
=======
from ml_peg.analysis.utils.utils import get_table_style
from math import floor

# Shared style rules to decorate Rank column like medals
RANK_STYLE_RULES = [
    {
        "if": {"filter_query": "{Rank} = 1", "column_id": "Rank"},
        "background": "linear-gradient(135deg, #fbbf24, #f59e0b)",
        "color": "#ffffff",
        "fontWeight": "700",
        "textAlign": "center",
        "borderRadius": "9999px",
    },
    {
        "if": {"filter_query": "{Rank} = 2", "column_id": "Rank"},
        "background": "linear-gradient(135deg, #d1d5db, #9ca3af)",
        "color": "#ffffff",
        "fontWeight": "700",
        "textAlign": "center",
        "borderRadius": "9999px",
    },
    {
        "if": {"filter_query": "{Rank} = 3", "column_id": "Rank"},
        "background": "linear-gradient(135deg, #fb923c, #f97316)",
        "color": "#ffffff",
        "fontWeight": "700",
        "textAlign": "center",
        "borderRadius": "9999px",
    },
]

def _rank_gradient_styles(data: list[dict]) -> list[dict]:
    """Generate background shades for Rank column from yellow→purple.

    Yellow (#fbbf24) for best (rank=1) to purple (#8b5cf6) for worst.
    """
    # Collect all numeric ranks present
    ranks = [row.get("Rank") for row in data if isinstance(row.get("Rank"), (int, float))]
    if not ranks:
        return []
    vmin, vmax = min(ranks), max(ranks)
    if vmax <= vmin:
        vmax = vmin + 1

    def lerp(a: float, b: float, t: float) -> float:
        return a + (b - a) * t

    def hex_to_rgb(h: str) -> tuple[int, int, int]:
        h = h.lstrip('#')
        return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)

    def rgb_to_hex(r: int, g: int, b: int) -> str:
        return f"#{r:02x}{g:02x}{b:02x}"

    # Yellow to purple
    y_rgb = hex_to_rgb("#fbbf24")
    p_rgb = hex_to_rgb("#8b5cf6")

    styles: list[dict] = []
    for val in ranks:
        t = (val - vmin) / (vmax - vmin)
        r = floor(lerp(y_rgb[0], p_rgb[0], t))
        g = floor(lerp(y_rgb[1], p_rgb[1], t))
        b = floor(lerp(y_rgb[2], p_rgb[2], t))
        color = rgb_to_hex(r, g, b)
        styles.append({
            "if": {"filter_query": f"{{Rank}} = {val}", "column_id": "Rank"},
            "backgroundColor": color,
            "color": "#000000",
        })
    return styles
>>>>>>> Stashed changes


def rebuild_table(filename: str | Path, id: str) -> DataTable:
    """
    Rebuild saved dash table.

    Parameters
    ----------
    filename
        Name of json file with saved table data.
    id
        ID for table.

    Returns
    -------
    DataTable
        Loaded Dash DataTable.

    Raises
    ------
    ValueError
        If the table JSON omits required ``thresholds`` metadata.
    """
    # Load JSON file
    with open(filename) as f:
        table_json = json.load(f)

    data = table_json["data"]
    columns = table_json["columns"]

    for column in columns:
        column_id = column.get("id")
        if column_id is None:
            continue
        if column_id == "Rank":
            column["type"] = "numeric"
            column.setdefault("format", rank_format())
        elif column.get("type") == "numeric" or is_numeric_column(data, column_id):
            column["type"] = "numeric"
            column.setdefault("format", sig_fig_format())
    tooltip_header = table_json["tooltip_header"]

    scored_data = calc_metric_scores(data)
    style = get_table_style(data, scored_data=scored_data)
    column_widths = calculate_column_widths([cols["name"] for cols in columns])

    style_cell_conditional: list[dict[str, object]] = []
    for column_id, width in column_widths.items():
        if width is None:
            continue
        col_width = f"{width}px"
        alignment = "left" if column_id == "MLIP" else "center"
        style_cell_conditional.append(
            {
                "if": {"column_id": column_id},
                "width": col_width,
                "minWidth": col_width,
                "maxWidth": col_width,
                "textAlign": alignment,
            }
        )

    table = DataTable(
        data=data,
        columns=columns,
        tooltip_header=tooltip_header,
        tooltip_delay=100,
        tooltip_duration=None,
        editable=True,
        id=id,
<<<<<<< Updated upstream
        style_data_conditional=style,
        style_cell_conditional=style_cell_conditional,
        sort_action="native",
        persistence=True,
        persistence_type="session",
        persisted_props=["data"],
=======
        style_data_conditional=style + _rank_gradient_styles(data) + RANK_STYLE_RULES,
        sort_action="native",
        style_table={"overflowX": "auto"},
        style_cell={
            "textAlign": "left",
            "padding": "1rem 1.5rem",
            "fontFamily": "var(--font-sans)",
            "fontSize": "0.875rem",
            "borderBottom": "1px solid #ffffff",
            "color": "#ffffff",
            "backgroundColor": "#000000",
        },
        style_cell_conditional=[
            {
                "if": {"column_id": "MLIP"},
                "color": "var(--color-primary)",
                "fontWeight": "600",
            },
            {
                "if": {"column_id": "Rank"},
                "textAlign": "center",
                "width": "80px",
                "minWidth": "80px",
                "maxWidth": "100px",
                "color": "#000000",
            },
        ],
        style_header={
            "backgroundColor": "#000000",
            "fontWeight": "600",
            "borderBottom": "2px solid #ffffff",
            "color": "#ffffff",
        },
        css=[
            {
                "selector": ".dash-table-container",
                "rule": "font-family: var(--font-sans);",
            }
        ],
>>>>>>> Stashed changes
    )

    thresholds = clean_thresholds(table_json.get("thresholds"))
    if thresholds is None or not thresholds:
        raise ValueError(f"No thresholds defined in table JSON: {filename}")

    table.thresholds = thresholds

    return table


def read_plot(filename: str | Path, id: str = "figure-1") -> Graph:
    """
    Read preprepared plotly Figure.

    Parameters
    ----------
    filename
        Name of json file with saved plot data.
    id
        ID for plot.

    Returns
    -------
    Graph
        Loaded plotly Graph.
    """
    fig = read_json(filename)
    # Apply a design-system-aligned theme
    fig.update_layout(
        font=dict(family="-apple-system, BlinkMacSystemFont, 'Inter', 'Segoe UI', Roboto, sans-serif", size=14, color="#0f172a"),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(gridcolor="#f1f5f9", linecolor="#e2e8f0", zerolinecolor="#e2e8f0"),
        yaxis=dict(gridcolor="#f1f5f9", linecolor="#e2e8f0", zerolinecolor="#e2e8f0"),
        hoverlabel=dict(bgcolor="#ffffff", font_size=12, font_family="-apple-system, BlinkMacSystemFont, 'Inter', 'Segoe UI', Roboto, sans-serif"),
    )
    return Graph(id=id, figure=fig)
