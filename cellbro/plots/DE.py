import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate
from plotly.subplots import make_subplots

import scout
# import cellbro.util.Components import Components
from cellbro.util.DashPage import DashPage
from cellbro.util.Param import *

figure_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=0, b=0, l=0, r=0),
)

de_params = ParamsDict(
    [
        Param(
            "method",
            "Method",
            default="t-test",
            type=str,
            description="",
            allowed_values=["logreg", "t-test", "wilcoxon", "t-test_overestim_var"],
        ),
        Param(
            "corr_method",
            "Correction Method",
            default="benjamini-hochberg",
            type=str,
            description="",
            allowed_values=[
                "benjamini-hochberg",
                "t-test",
                "t-test_overestim_var",
                "wilcoxon",
            ],
        ),
    ]
)


def apply(dataset, params):
    if params["submit"] is None:
        groupby_options = get_groupby_options(dataset=dataset)
        if len(groupby_options) == 0:
            return [[], None, [], None, {"display": "none"}, {"display": "none"}]

        groupby_selected = groupby_options[0]
        refs_options = get_reference_options(dataset=dataset, groupby=groupby_selected)
        refs_selected = next(
            iter(dataset.adata.uns[f"rank_genes_{groupby_selected}"].keys())
        )
        return [
            groupby_options,
            groupby_selected,
            refs_options,
            refs_selected,
            {"display": "block"},
            {"display": "block"},
        ]

    key = f"rank_genes_{params['groupby']}"

    if key not in dataset.adata.uns.keys():
        scout.tl.rank_marker_genes(dataset.adata, groupby=params["groupby"])

    groupby_options = get_groupby_options(dataset=dataset)
    refs_options = get_reference_options(dataset, params["groupby"])

    groupby_selected = params["groupby"]
    refs_selected = next(iter(dataset.adata.uns[key].keys()))
    return [
        groupby_options,
        groupby_selected,
        refs_options,
        refs_selected,
        {"display": "block"},
        {"display": "block"},
    ]


def plot_pval_histogram(dataset, params):
    key = f"rank_genes_{params['groupby']}"
    fig = scout.ply.pval_histogram(
        dataset.adata.uns[key][params["reference"]], layout=figure_layout
    )
    return fig


def plot_de_volcano(dataset, params):
    key = f"rank_genes_{params['groupby']}"
    if params["reference"] is None:
        params["reference"] = list(dataset.adata.uns[key].keys())[0]
    fig = scout.ply.marker_volcano(
        dataset.adata.uns[key][params["reference"]], layout=figure_layout
    )
    return fig


def get_groupby_options(dataset):
    return [x[11:] for x in dataset.adata.uns.keys() if x.startswith("rank_genes_")]


def get_reference_options(dataset, groupby):
    groupby = dataset.adata.uns.get(f"rank_genes_{groupby}", None)
    if groupby is None:
        return []
    return dict([(key, f"{key} vs. Rest") for key in groupby.keys()])
