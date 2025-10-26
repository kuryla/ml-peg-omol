"""Fixtures for MLIP results analysis."""

from __future__ import annotations

from collections.abc import Callable
import functools
from json import dump
from pathlib import Path
from typing import Any

from dash import dash_table
import numpy as np
import plotly.graph_objects as go

from ml_peg.analysis.utils.utils import calc_ranks, calc_scores


def plot_parity(
    title: str | None = None,
    x_label: str | None = None,
    y_label: str | None = None,
    hoverdata: dict | None = None,
    filename: str = "parity.json",
) -> Callable:
    """
    Plot parity plot of MLIP results against reference data.

    Parameters
    ----------
    title
        Graph title.
    x_label
        Label for x-axis. Default is `None`.
    y_label
        Label for y-axis. Default is `None`.
    hoverdata
        Hover data dictionary. Default is `{}`.
    filename
        Filename to save plot as JSON. Default is "parity.json".

    Returns
    -------
    Callable
        Decorator to wrap function.
    """

    def plot_parity_decorator(func: Callable) -> Callable:
        """
        Decorate function to plot parity.

        Parameters
        ----------
        func
            Function being wrapped.

        Returns
        -------
        Callable
            Wrapped function.
        """

        @functools.wraps(func)
        def plot_parity_wrapper(*args, **kwargs) -> dict[str, Any]:
            """
            Wrap function to plot parity.

            Parameters
            ----------
            *args
                Arguments to pass to the function being wrapped.
            **kwargs
                Key word arguments to pass to the function being wrapped.

            Returns
            -------
            dict
                Results dictionary.
            """
            results = func(*args, **kwargs)
            ref = results["ref"]

            hovertemplate = "<b>Pred: </b>%{x}<br>" + "<b>Ref: </b>%{y}<br>"

            customdata = []
            if hoverdata:
                for i, key in enumerate(hoverdata):
                    hovertemplate += f"<b>{key}: </b>%{{customdata[{i}]}}<br>"
                customdata = list(zip(*hoverdata.values(), strict=True))

            fig = go.Figure()
            for mlip, value in results.items():
                if mlip == "ref":
                    continue
                fig.add_trace(
                    go.Scatter(
                        x=value,
                        y=ref,
                        name=mlip,
                        mode="markers",
                        customdata=customdata,
                        hovertemplate=hovertemplate,
                    )
                )

            # Compute plot limits directly from data (no Kaleido dependency)
            x_values: list[float] = []
            for mlip, value in results.items():
                if mlip == "ref":
                    continue
                x_values.extend(value)

            if x_values and ref:
                x_min = min(min(x_values), min(ref))
                x_max = max(max(x_values), max(ref))
            else:
                # Fallback if no data
                x_min, x_max = 0.0, 1.0

            lims = [x_min, x_max]

            fig.add_trace(
                go.Scatter(
                    x=lims,
                    y=lims,
                    mode="lines",
                    showlegend=False,
                )
            )

            fig.update_layout(
                title={"text": title},
                xaxis={"title": {"text": x_label}},
                yaxis={"title": {"text": y_label}},
            )

            fig.update_traces()

            # Write to file
            Path(filename).parent.mkdir(parents=True, exist_ok=True)
            fig.write_json(filename)

            return results

        return plot_parity_wrapper

    return plot_parity_decorator


def plot_scatter(
    title: str | None = None,
    x_label: str | None = None,
    y_label: str | None = None,
    show_line: bool = False,
    hoverdata: dict | None = None,
    filename: str = "scatter.json",
) -> Callable:
    """
    Plot scatter plot of MLIP results.

    Parameters
    ----------
    title
        Graph title.
    x_label
        Label for x-axis. Default is `None`.
    y_label
        Label for y-axis. Default is `None`.
    show_line
        Whether to show line between points. Default is False.
    hoverdata
        Hover data dictionary. Default is `{}`.
    filename
        Filename to save plot as JSON. Default is "scatter.json".

    Returns
    -------
    Callable
        Decorator to wrap function.
    """

    def plot_scatter_decorator(func: Callable) -> Callable:
        """
        Decorate function to plot scatter.

        Parameters
        ----------
        func
            Function being wrapped.

        Returns
        -------
        Callable
            Wrapped function.
        """

        @functools.wraps(func)
        def plot_scatter_wrapper(*args, **kwargs) -> dict[str, Any]:
            """
            Wrap function to plot scatter.

            Parameters
            ----------
            *args
                Arguments to pass to the function being wrapped.
            **kwargs
                Key word arguments to pass to the function being wrapped.

            Returns
            -------
            dict
                Results dictionary.
            """
            results = func(*args, **kwargs)

            hovertemplate = "<b>Pred: </b>%{x}<br>" + "<b>Ref: </b>%{y}<br>"
            customdata = []
            if hoverdata:
                for i, key in enumerate(hoverdata):
                    hovertemplate += f"<b>{key}: </b>%{{customdata[{i}]}}<br>"
                customdata = list(zip(*hoverdata.values(), strict=True))

            mode = "lines+markers" if show_line else "markers"

            fig = go.Figure()
            for mlip, value in results.items():
                name = "Reference" if mlip == "ref" else mlip
                fig.add_trace(
                    go.Scatter(
                        x=value[0],
                        y=value[1],
                        name=name,
                        mode=mode,
                        customdata=customdata,
                        hovertemplate=hovertemplate,
                    )
                )

            fig.update_layout(
                title={"text": title},
                xaxis={"title": {"text": x_label}},
                yaxis={"title": {"text": y_label}},
            )

            fig.update_traces()

            # Write to file
            Path(filename).parent.mkdir(parents=True, exist_ok=True)
            fig.write_json(filename)

            return results

        return plot_scatter_wrapper

    return plot_scatter_decorator


def build_table(
    filename: str = "table.json",
    metric_tooltips: dict[str, str] | None = None,
) -> Callable:
    """
    Build table MLIP results.

    Parameters
    ----------
    filename
        Filename to save table.
    metric_tooltips
        Tooltips for table metric headers.

    Returns
    -------
    Callable
        Decorator to wrap function.
    """

    def build_table_decorator(func: Callable) -> Callable:
        """
        Decorate function to build table.

        Parameters
        ----------
        func
            Function being wrapped.

        Returns
        -------
        Callable
            Wrapped function.
        """

        @functools.wraps(func)
        def build_table_wrapper(*args, **kwargs) -> dict[str, Any]:
            """
            Wrap function to build table.

            Parameters
            ----------
            *args
                Arguments to pass to the function being wrapped.
            **kwargs
                Key word arguments to pass to the function being wrapped.

            Returns
            -------
            dict
                Results dictionary.
            """
            results = func(*args, **kwargs)
            # Form of results is
            # results = {
            #     metric_1: {mlip_1: value_1, mlip_2: value_2},
            #     metric_2: {mlip_1: value_3, mlip_2: value_4},
            # }

            metrics_columns = ("MLIP",) + tuple(results.keys())
            # Use MLIP keys from first (any) metric keys
            mlips = next(iter(results.values())).keys()

            metrics_data = []
            for mlip in mlips:
                metrics_data.append(
                    {"MLIP": mlip}
                    | {key: value[mlip] for key, value in results.items()}
                    | {"id": mlip},
                )

            summary_tooltips = {
                "MLIP": "Name of the model",
                "Score": "Average of metrics (lower is better)",
                "Rank": "Model rank based on score (lower is better)",
            }
            if metric_tooltips:
                tooltip_header = metric_tooltips | summary_tooltips
            else:
                tooltip_header = summary_tooltips

            metrics_data = calc_scores(metrics_data)
            metrics_data = calc_ranks(metrics_data)
            metrics_columns += ("Score", "Rank")

            table = dash_table.DataTable(
                metrics_data,
                [{"name": i, "id": i} for i in metrics_columns if i != "id"],
                id="metrics",
                tooltip_header=tooltip_header,
            )

            # Save dict of table to be loaded
            Path(filename).parent.mkdir(parents=True, exist_ok=True)
            with open(filename, "w") as fp:
                dump(
                    {
                        "data": table.data,
                        "columns": table.columns,
                        "tooltip_header": tooltip_header,
                    },
                    fp,
                )

            return results

        return build_table_wrapper

    return build_table_decorator
