import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import Input, Output, State, dcc, html

from cellbro.util.DashFigure import DashFigure
from cellbro.util.DashAction import DashAction
import cellbro.util.Components as Components
from cellbro.util.Param import *

import scout

violin_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=5, b=5, l=5, r=5),
)

class PlotViolin(DashAction):
    def plot(self, feature, groupby):
        fig = scout.ply.violin(self.dataset.adata, y=feature, groupby=groupby, layout=violin_layout)
        return [fig]

    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-violin-plot", component_property="figure"),
        ]

        inputs = {
            "feature": Input(component_id=f"{self.page_id_prefix}-violin-feature", component_property="value"),
            "groupby": Input(component_id=f"{self.page_id_prefix}-violin-groupby", component_property="value"),
        }

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(feature, groupby):
            return self.plot(feature, groupby)


class Violin(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.actions.update(
            plot_violin=PlotViolin(dataset, self.page_id_prefix)
        )

    def create_layout(self) -> list:
        var_names = [(f, f) for f in sorted(list(self.dataset.adata.var_names))]

        other = [
            (f, f.replace("_", " ").capitalize())
            for f in sorted(list(self.dataset.get_numeric()))
        ]

        features = dict(other + var_names)

        groupbys = dict(
            [
                (k, k.replace("_", " ").capitalize())
                for k in self.dataset.get_categoric()
            ]
        )

        type_params = Components.FigureParamTab(self.page_id_prefix, tab_label="Type", children=[
            # Features
            html.Div([
                html.Label("Feature"),
                dcc.Dropdown(
                    features,
                    value=list(features.keys())[0],
                    id=f"{self.page_id_prefix}-violin-feature",
                    clearable=False,
                )
            ], className="param-row-stacked"),
            # Groupby select
            html.Div([
                html.Label("Group By"),
                dcc.Dropdown(
                    groupbys,
                    value=None,
                    id=f"{self.page_id_prefix}-violin-groupby",
                    clearable=True,
                ),
            ], className="param-row-stacked")
        ])

        figure_params = Components.FigureParams(self.page_id_prefix, tabs=[type_params])


        figure = html.Div(
            children=[
                html.Div(
                    children=figure_params.create_layout(), className="fig-header"
                ),
                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id=f"{self.page_id_prefix}-violin-plot", className="secondary-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    className="secondary-body",
                ),
            ],
            className="secondary",
        )

        return figure

    def get_sidebar_params(self) -> list:
        return []
