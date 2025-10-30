"""Build main Dash application."""

from __future__ import annotations

from importlib import import_module
import warnings

from dash import Dash, Input, Output, callback, html
from dash.dash_table import DataTable
from dash.dcc import Store, Tab, Tabs, Location
from dash.html import H1, H3, Div
from yaml import safe_load
from ml_peg.models import MODELS_ROOT

<<<<<<< Updated upstream
from ml_peg.analysis.utils.utils import calc_ranks, calc_table_scores, get_table_style
=======
from ml_peg.analysis.utils.utils import calc_ranks, calc_scores, get_table_style
from ml_peg.app.utils.load import RANK_STYLE_RULES
>>>>>>> Stashed changes
from ml_peg.app import APP_ROOT
from ml_peg.app.utils.build_components import build_weight_components
from ml_peg.app.utils.register_callbacks import register_benchmark_to_category_callback
from ml_peg.app.utils.utils import calculate_column_widths, rank_format, sig_fig_format
from ml_peg.models.get_models import get_model_names
from ml_peg.models.models import current_models

MODELS = get_model_names(current_models)


def get_all_tests(
    category: str = "*",
) -> tuple[
    dict[str, dict[str, list[Div]]],
    dict[str, dict[str, DataTable]],
    dict[str, dict[str, dict]],
]:
    """
    Get layout and register callbacks for all categories.

    Parameters
    ----------
    category
        Name of category directory to search for tests. Default is '*'.

    Returns
    -------
    tuple[dict[str, dict[str, list[Div]]], dict[str, dict[str, DataTable]]]
        Layouts and tables for all categories.
    """
    # Find Python files e.g. app_OC157.py in mlip_tesing.app module.
    # We will get the category from the parent's parent directory
    # E.g. ml_peg/app/surfaces/OC157/app_OC157.py -> surfaces
    tests = APP_ROOT.glob(f"{category}/*/app*.py")
    layouts = {}
    tables = {}
    infos: dict[str, dict[str, dict]] = {}

    # Build all layouts, and register all callbacks to main app.
    for test in tests:
        try:
            # Import test layout/callbacks
            test_name = test.parent.name
            category_name = test.parent.parent.name
            test_module = import_module(
                f"ml_peg.app.{category_name}.{test_name}.app_{test_name}"
            )
            test_app = test_module.get_app()

            # Get layouts and tables for each category/test
            if category_name not in layouts:
                layouts[category_name] = {}
                tables[category_name] = {}
                infos[category_name] = {}
            layouts[category_name][test_app.name] = test_app.layout
            tables[category_name][test_app.name] = test_app.table
            infos[category_name][test_app.name] = {
                "description": getattr(test_app, "description", ""),
                "docs_url": getattr(test_app, "docs_url", None),
            }
        except FileNotFoundError as err:
            warnings.warn(
                f"Unable to load layout for {test_name} in {category_name} category. "
                f"Full error:\n{err}",
                stacklevel=2,
            )
            continue

        # Register test callbacks
        try:
            test_app.register_callbacks()
        except FileNotFoundError as err:
            warnings.warn(
                f"Unable to register callbacks for {test_name} in {category_name} "
                f"category. Full error:\n{err}",
                stacklevel=2,
            )
            continue

    return layouts, tables, infos


def build_category(
    all_layouts: dict[str, dict[str, list[Div]]],
    all_tables: dict[str, dict[str, DataTable]],
) -> tuple[dict[str, list[Div]], dict[str, DataTable]]:
    """
    Build category layouts and summary tables.

    Parameters
    ----------
    all_layouts
        Layouts of all tests, grouped by category.
    all_tables
        Tables for all tests, grouped by category.

    Returns
    -------
    tuple[dict[str, list[Div]], dict[str, DataTable]]
        Dictionary of category layouts, and dictionary of category summary tables.
    """
    # Take all tables in category, build new table, and set layout
    category_layouts = {}
    category_tables = {}

    # `category` corresponds to the category's directory name
    # We will use the loaded `category_title` for IDs/dictionary keys returned
    for category in all_layouts:
        # Get category name and description
        try:
            with open(APP_ROOT / category / f"{category}.yml") as file:
                category_info = safe_load(file)
                category_title = category_info.get("title", category)
                category_descrip = category_info.get("description", "")
        except FileNotFoundError:
            category_title = category
            category_descrip = ""

        # Build category summary table
        summary_table = build_summary_table(
            all_tables[category], table_id=f"{category_title}-summary-table"
        )
        category_tables[category_title] = summary_table

        # Build weight components for category summary table
        weight_components = build_weight_components(
            header="Benchmark weights",
            table=summary_table,
            column_widths=getattr(summary_table, "column_widths", None),
        )

        # Build full layout with summary table, weight controls, and test layouts
        category_layouts[category_title] = Div(
            [
                H1(category_title),
                H3(category_descrip),
                summary_table,
                Store(
                    id=f"{category_title}-summary-table-computed-store",
                    storage_type="session",
                    data=summary_table.data,
                ),
                weight_components,
                Div([all_layouts[category][test] for test in all_layouts[category]]),
            ]
        )

        # Register benchmark table -> category table callbacks
        for test_name, benchmark_table in all_tables[category].items():
            register_benchmark_to_category_callback(
                benchmark_table_id=benchmark_table.id,
                category_table_id=f"{category_title}-summary-table",
                benchmark_column=test_name,
            )

    return category_layouts, category_tables


def build_summary_table(
    tables: dict[str, DataTable], table_id: str = "summary-table"
) -> DataTable:
    """
    Build summary table from a set of tables.

    Parameters
    ----------
    tables
        Dictionary of tables to be summarised.
    table_id
        ID of table being built. Default is 'summary-table'.

    Returns
    -------
    DataTable
        Summary table with score from table being summarised.
    """
    summary_data = {}
    for category_name, table in tables.items():
        # Prepare rows for all current models
        if not summary_data:
            summary_data = {model: {} for model in MODELS}

        for row in table.data:
            # Category tables may include models not to be included
            if row["MLIP"] in summary_data:
                summary_data[row["MLIP"]][category_name] = row["Score"]

    data = []
    for mlip in summary_data:
        data.append({"MLIP": mlip} | summary_data[mlip])

    data = calc_table_scores(data)
    data = calc_ranks(data)

    columns_headers = ("MLIP",) + tuple(tables.keys()) + ("Score", "Rank")

    columns = [{"name": headers, "id": headers} for headers in columns_headers]
    for column in columns:
        column_id = column["id"]
        if column_id == "Rank":
            column["type"] = "numeric"
            column["format"] = rank_format()
        elif column_id != "MLIP":
            column["type"] = "numeric"
            column["format"] = sig_fig_format()

    style = get_table_style(data)

    # Calculate column widths based on column names
    column_widths = calculate_column_widths(columns_headers)
    style_cell_conditional = []
    for column_id, width in column_widths.items():
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
        id=table_id,
        sort_action="native",
<<<<<<< Updated upstream
        style_data_conditional=style,
        style_cell_conditional=style_cell_conditional,
        persistence=True,
        persistence_type="session",
        persisted_props=["data"],
=======
        style_data_conditional=style + RANK_STYLE_RULES,
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
    table.column_widths = column_widths
    return table


def build_tabs(
    full_app: Dash,
    layouts: dict[str, list[Div]],
    summary_table: DataTable,
    weight_components: Div,
) -> None:
    """
    Build tab layouts and summary tab.

    Parameters
    ----------
    full_app
        Full application with all sub-apps.
    layouts
        Layouts for all tabs.
    summary_table
        Summary table with score from each category.
    weight_components
        Weight sliders, text boxes and reset button.
    """
    all_tabs = [Tab(label="Summary", value="summary-tab", id="summary-tab")] + [
        Tab(label=category_name, value=category_name) for category_name in layouts
    ]

    # Header (site-wide)
    header = html.Header(
        className="header",
        children=[
            html.Div(
                className="header-content",
                children=[
                    html.A(className="logo", href="/", children=[html.Span("ML-PEG", className="logo-text")]),
                    html.Nav(
                        className="nav",
                        id="site-nav",
                        children=[
                            html.A("Leaderboard", href="/", className="nav-link active"),
                            html.A("Models", href="/models", className="nav-link"),
                            html.A("Datasets", href="/datasets", className="nav-link"),
                            html.A("Submit Model", href="#", className="button button-primary"),
                        ],
                    ),
                ],
            )
        ],
    )

    # Hero section
    # No hero on dark minimalist layout

    # Home/Leaderboard content creator
    def leaderboard_content() -> html.Div:
        return html.Div(
            className="container",
            children=[
                html.Div(
                    className="card card--borderless",
                    children=[
                        html.Div(
                            className="card-header",
                            children=[html.H3("MLIP Benchmarking", className="card-title")],
                        ),
                        html.Div(
                            className="card-body",
                            children=[
                                Tabs(
                                    id="all-tabs",
                                    value="summary-tab",
                                    children=all_tabs,
                                    className="dash-tabs",
                                ),
                                Div(id="tabs-content"),
                            ],
                        ),
                    ],
                )
            ],
        )

    # Footer
    footer = html.Footer(
        className="footer",
        children=[
            html.Div(
                className="container",
                children=[
                    html.Div(
                        className="footer-content",
                        children=[
                            html.Div(
                                className="footer-section",
                                children=[
                                    html.H4("ML-PEG"),
                                    html.P(
                                        "A comprehensive benchmarking framework for machine learning interatomic potentials."
                                    ),
                                ],
                            ),
                            html.Div(
                                className="footer-section",
                                children=[
                                    html.H4("Resources"),
                                    html.Ul(
                                        className="footer-links",
                                        children=[
                                            html.Li(html.A("Documentation", href="#")),
                                            html.Li(html.A("GitHub Repository", href="#")),
                                            html.Li(html.A("Contributing Guide", href="#")),
                                        ],
                                    ),
                                ],
                            ),
                            html.Div(
                                className="footer-section",
                                children=[
                                    html.H4("Community"),
                                    html.Ul(
                                        className="footer-links",
                                        children=[
                                            html.Li(html.A("Discussion Forum", href="#")),
                                            html.Li(html.A("Issue Tracker", href="#")),
                                            html.Li(html.A("Twitter", href="#")),
                                        ],
                                    ),
                                ],
                            ),
                            html.Div(
                                className="footer-section",
                                children=[
                                    html.H4("Legal"),
                                    html.Ul(
                                        className="footer-links",
                                        children=[
                                            html.Li(html.A("Terms of Use", href="#")),
                                            html.Li(html.A("Privacy Policy", href="#")),
                                            html.Li(html.A("License", href="#")),
                                        ],
                                    ),
                                ],
                            ),
                        ],
                    ),
                    html.Div(
                        className="footer-bottom",
                        children=[
                            html.P("© 2025 ML-PEG. All rights reserved. | Hosted by STFC"),
                        ],
                    ),
                ],
            )
        ],
    )

    # Router: page content goes here based on URL
    full_app.layout = html.Div(
        [
            header,
            Location(id="url", refresh=False),
            html.Div(id="page-root", children=[leaderboard_content()]),
            footer,
        ]
    )

    @callback(Output("tabs-content", "children"), Input("all-tabs", "value"))
    def select_tab(tab) -> Div:
        """
        Select tab contents to be displayed.

        Parameters
        ----------
        tab
            Name of tab selected.

        Returns
        -------
        Div
            Summary or tab contents to be displayed.
        """
        if tab == "summary-tab":
            return Div(
                className="card",
                children=[
                    html.Div(className="card-header", children=[html.H3("Benchmarks Summary", className="card-title")]),
                    html.Div(className="card-body", children=[summary_table, weight_components]),
                    Store(id="summary-table-scores-store", storage_type="session"),
                ],
            )

        # CATEGORY PAGE: simple card with the original layout
        return Div(className="card", children=[
            html.Div(className="card-header", children=[html.H3(tab, className="card-title")]),
            html.Div(className="card-body", children=[layouts[tab]]),
        ])

    # Router: swap page content by URL
    @callback(Output("page-root", "children"), Input("url", "pathname"))
    def render_page(pathname):
        if pathname in ("/", "", None):
            # Home: leaderboard tabs only
            return [leaderboard_content()]
        if pathname == "/models":
            # Models page
            import yaml
            with open(MODELS_ROOT / "models.yml") as f:
                all_models = yaml.safe_load(f) or {}
            rows = []
            for name, cfg in all_models.items():
                params = cfg.get("kwargs", {}) or {}
                param_str = ", ".join(f"{k}={v}" for k, v in params.items()) if params else ""
                rows.append({
                    "Name": f"[{name}](/models/{name})",
                    "Class": cfg.get("class_name", ""),
                    "Module": cfg.get("module", ""),
                    "Device": cfg.get("device", ""),
                    "Parameters": param_str,
                })
            model_table = DataTable(
                data=rows,
                columns=[
                    {"name": "Name", "id": "Name", "presentation": "markdown"},
                    {"name": "Class", "id": "Class"},
                    {"name": "Module", "id": "Module"},
                    {"name": "Device", "id": "Device"},
                    {"name": "Parameters", "id": "Parameters"},
                ],
                id="models-table",
                sort_action="native",
                style_table={"overflowX": "auto"},
                style_cell={
                    "textAlign": "left",
                    "padding": "1rem 1.5rem",
                    "fontFamily": "var(--font-sans)",
                    "fontSize": "0.875rem",
                    "borderBottom": "1px solid #ffffff",
                    "color": "var(--color-primary)",
                    "backgroundColor": "#000000",
                },
                style_header={
                    "backgroundColor": "#000000",
                    "fontWeight": "600",
                    "borderBottom": "2px solid #ffffff",
                    "color": "#ffffff",
                },
            )
            return html.Div(className="container", children=[
                html.Div(className="card card--borderless", children=[
                    html.Div(className="card-header", children=[html.H3("Models", className="card-title")]),
                    html.Div(className="card-body", children=[model_table]),
                ])
            ])
        if pathname and pathname.startswith("/models/"):
            # Model detail page
            import yaml
            name = pathname.split("/", 2)[-1]
            with open(MODELS_ROOT / "models.yml") as f:
                all_models = yaml.safe_load(f) or {}
            cfg = all_models.get(name)
            if not cfg:
                return html.Div(className="container", children=[
                    html.Div(className="card", children=[
                        html.Div(className="card-header", children=[html.H3("Model not found", className="card-title")]),
                        html.Div(className="card-body", children=[html.P(f"No model named '{name}' in models.yml.")]),
                    ])
                ])
            meta_items = [
                ("Name", name),
                ("Class", cfg.get("class_name", "")),
                ("Module", cfg.get("module", "")),
                ("Device", cfg.get("device", "")),
            ]
            params = cfg.get("kwargs", {}) or {}
            return html.Div(className="container", children=[
                html.Div(className="card", children=[
                    html.Div(className="card-header", children=[html.H3(f"Model: {name}", className="card-title")]),
                    html.Div(className="card-body", children=[
                        html.P("Configuration"),
                        html.Ul(children=[html.Li(html.Span(f"{k}: {v}")) for k, v in meta_items]),
                        html.Br(),
                        html.P("Parameters"),
                        html.Ul(children=[html.Li(html.Span(f"{k}: {v}")) for k, v in params.items()]) if params else html.P("No extra parameters"),
                        html.Br(),
                        html.A("← Back to Models", href="/models", className="nav-link"),
                    ]),
                ])
            ])
        if pathname == "/datasets":
            sections = []
            for category_dir in APP_ROOT.iterdir():
                if not category_dir.is_dir():
                    continue
                info_file = category_dir / f"{category_dir.name}.yml"
                if not info_file.exists():
                    continue
                try:
                    with open(info_file) as f:
                        info = safe_load(f) or {}
                except Exception:
                    info = {}
                title = info.get("title", category_dir.name)
                desc = info.get("description", "")
                datasets = []
                for sub in category_dir.iterdir():
                    if sub.is_dir() and (sub / f"app_{sub.name}.py").exists():
                        datasets.append(sub.name)
                datasets.sort()
                sections.append(
                    html.Div(className="card", children=[
                        html.Div(className="card-header", children=[html.H3(title, className="card-title")]),
                        html.Div(className="card-body", children=[
                            html.P(desc),
                            html.Div(
                                children=[html.A(ds, href=f"/datasets/{category_dir.name}/{ds}", className="badge badge-info") for ds in datasets],
                                style={"display": "flex", "gap": "10px", "flexWrap": "wrap"},
                            ),
                        ]),
                    ])
                )
            return html.Div(className="container", children=sections)
        if pathname and pathname.startswith("/datasets/"):
            parts = pathname.strip("/").split("/")
            if len(parts) >= 3:
                category, dataset = parts[1], parts[2]
                # Attempt to import module to read DOCS_URL
                docs_url = None
                try:
                    mod = __import__(f"ml_peg.app.{category}.{dataset}.app_{dataset}", fromlist=["*"])
                    docs_url = getattr(mod, "DOCS_URL", None)
                except Exception:
                    docs_url = None
                # Build data availability info
                from pathlib import Path
                data_dir = APP_ROOT / "data" / category / dataset
                available = data_dir.exists()
                metrics_files = list(data_dir.glob("*_metrics_table.json")) if available else []
                models_dirs = [p.name for p in data_dir.iterdir() if p.is_dir()] if available else []
                body_children = [
                    html.P(f"Category: {category}"),
                    html.P(f"Dataset: {dataset}"),
                ]
                if docs_url:
                    body_children.append(html.P(["Docs: ", html.A(docs_url, href=docs_url, target="_blank")]))
                body_children.append(html.P(f"Data directory: {'available' if available else 'missing'}"))
                if metrics_files:
                    body_children.append(html.P(f"Metrics files: {', '.join(p.name for p in metrics_files)}"))
                if models_dirs:
                    body_children.append(html.P("Models with structures:"))
                    body_children.append(html.Div(children=[html.Span(m, className="badge badge-primary") for m in sorted(models_dirs)], style={"display":"flex","gap":"8px","flexWrap":"wrap"}))
                body_children.append(html.Br())
                body_children.append(html.A("← Back to Datasets", href="/datasets", className="nav-link"))
                return html.Div(className="container", children=[
                    html.Div(className="card", children=[
                        html.Div(className="card-header", children=[html.H3(f"Dataset: {dataset}", className="card-title")]),
                        html.Div(className="card-body", children=body_children),
                    ])
                ])
            # fallback
            return html.Div(className="container", children=[html.Div(className="card", children=[
                html.Div(className="card-header", children=[html.H3("Dataset not found", className="card-title")]),
                html.Div(className="card-body", children=[html.P("The requested dataset page was not found.")]),
            ])])
        # 404
        return html.Div(className="container", children=[html.Div(className="card", children=[
            html.Div(className="card-header", children=[html.H3("Page not found", className="card-title")]),
            html.Div(className="card-body", children=[html.P("The requested page was not found.")]),
        ])])


def build_full_app(full_app: Dash, category: str = "*") -> None:
    """
    Build full app layout and register callbacks.

    Parameters
    ----------
    full_app
        Full application with all sub-apps.
    category
        Category to build app for. Default is `*`, corresponding to all categories.
    """
    # Get layouts and tables for each test, grouped by categories
    all_layouts, all_tables, all_infos = get_all_tests(category=category)

    if not all_layouts:
        raise ValueError("No tests were built successfully")

    # Combine tests into categories and create category summary
    category_layouts, category_tables = build_category(all_layouts, all_tables)
    # Build overall summary table
    summary_table = build_summary_table(category_tables)
    weight_components = build_weight_components(
        header="Category weights",
        table=summary_table,
        column_widths=getattr(summary_table, "column_widths", None),
    )
    # Build summary and category tabs
    build_tabs(full_app, category_layouts, summary_table, weight_components)
