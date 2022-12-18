import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import Input, Output, State, dcc, html

from cellbro.util.Param import *

import scout

violin_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=5, b=5, l=5, r=5),
)


class Violin:
    def __init__(self, dataset):
        self.dataset = dataset

    def plot(self, feature, groupby):
        fig = scout.ply.violin(self.dataset.adata, y=feature, groupby=groupby, layout=violin_layout)

        return [fig]

    @staticmethod
    def create_layout(dataset):
        var_names = [(f, f) for f in sorted(list(dataset.adata.var_names))]

        other = [
            (f, f.replace("_", " ").capitalize())
            for f in sorted(list(dataset.get_numeric()))
        ]

        features = dict(other + var_names)

        groupbys = dict(
            [
                (k, k.replace("_", " ").capitalize())
                for k in dataset.get_categoric()
            ]
        )
        figure = html.Div(
            children=[
                html.Div(
                    children=[
                        # Features
                        html.Div(
                            children=[
                                html.Label("Feature"),
                                dcc.Dropdown(
                                    features,
                                    value=list(features.keys())[0],
                                    id="violin-feature",
                                    clearable=False,
                                ),
                            ],
                            className="param-column",
                        ),
                        # Groupby select
                        html.Div(
                            children=[
                                html.Label("Group By"),
                                dcc.Dropdown(
                                    groupbys,
                                    value=None,
                                    id="violin-groupby",
                                    clearable=True,
                                ),
                            ],
                            className="param-column",
                        ),
                    ],
                    id="violin-select",
                    className="secondary-select top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            id="violin-projection",
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id="violin-plot", className="secondary-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    id="violin-figure",
                    className="secondary-figure",
                ),
            ],
            className="secondary",
        )

        return figure
