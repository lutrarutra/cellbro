import plotly.graph_objects as go
import plotly.express as px
from dash import html, dcc, Input, Output, State
import dash_bootstrap_components as dbc
from plotly.subplots import make_subplots

from cellbro.util.Param import *

import scanpy as sc
import pandas as pd
import numpy as np

figure_layout = go.Layout(
    paper_bgcolor='white',
    plot_bgcolor='white',
    xaxis=dict(showgrid=False, zeroline=False, visible=False, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=0, b=0, l=0, r=0),
)

class QC():
    def __init__(self, dataset):
        self.dataset = dataset

    def plot(self):
        fig = make_subplots(rows=1, cols=3)

        violins = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
        x = self.dataset.adata.obs[violins]

        df = pd.DataFrame(
            x, columns=violins
        )

        for i, feature in enumerate(violins):
            fig.add_trace(
                    go.Violin(
                        y=df[feature],
                        line_color=sc.pl.palettes.default_20[i],
                        box_visible=True,
                        points="all",
                        pointpos=0,
                        marker=dict(size=2),
                        jitter=0.6,
                ), row=1, col=i+1
            )
            fig["layout"][f"xaxis{i+1 if i != 0 else ''}"]["title"] = feature.replace("_", " ").title()

        fig.update_layout(figure_layout)

        return [fig]
    
        
    @staticmethod
    def get_callback_outputs():
        return [
            Output(component_id="bottom-plot", component_property="figure"),
        ]
        
    # Inputs to Projection
    @staticmethod
    def get_callback_inputs():
        return {
            "submit": Input(component_id="qc-submit", component_property="n_clicks"),
            # "feature": Input(component_id="qc-violin-feature", component_property="value"),
            # "groupby": Input(component_id="violin-groupby", component_property="value"),
        }

    @staticmethod
    def create_layout(dataset):
        sidebar = html.Div(children=[
            dcc.Loading(type="circle", children=[
                html.Div([
                    dbc.Button("Plot", color="primary", className="mr-1", id="qc-submit"),
                ], id="bottom-sidebar-footer")
            ],)
        ], id="bottom-sidebar")

        figure = html.Div(children=[
            dcc.Loading(
                id="loading-bottom", type="circle",
                children=[html.Div(dcc.Graph(id="bottom-plot"))],
            )
        ], id="bottom-figure")

        return sidebar, figure