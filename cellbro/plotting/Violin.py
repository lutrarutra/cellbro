import plotly.graph_objects as go
import plotly.express as px
from dash import html, dcc, Input, Output, State
import dash_bootstrap_components as dbc

from cellbro.util.Param import *

import scanpy as sc
import pandas as pd
import numpy as np

violin_layout = go.Layout(
    paper_bgcolor='white',
    plot_bgcolor='white',
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=0, b=0, l=0, r=0),
)

class Violin():
    def __init__(self, dataset):
        self.dataset = dataset

    def plot(self, feature, groupby):
        fig = go.Figure()

        if feature not in self.dataset.adata.obs_keys():
            x = self.dataset.adata[:, feature].X.toarray().T.squeeze()
        else:
            x = self.dataset.adata.obs[feature].values.squeeze()

        if groupby is not None:
            f = self.dataset.adata.obs[groupby].values
            df = pd.DataFrame(
                np.array([x,f]).T, columns=[feature, groupby]
            )
            groups = df[groupby].unique()
            for i, g in enumerate(groups):
                fig.add_trace(go.Violin(
                    x=df[df[groupby] == g][groupby],
                    y=df[df[groupby] == g][feature],
                    name=g,
                    line_color=sc.pl.palettes.default_20[i],
                    box_visible=True,
                    points="all",
                    pointpos=0,
                    marker=dict(size=2),
                    jitter=0.6,
                ))
            violin_layout["legend"] = dict(title=groupby.capitalize())
            violin_layout["xaxis"]["showticklabels"] = True
            violin_layout["xaxis_title"] = groupby.capitalize()
            
        else:
            df = pd.DataFrame(
                x, columns=[feature]
            )
            fig.add_trace(go.Violin(
                y=df[feature],
                line_color=sc.pl.palettes.default_20[0],
                box_visible=True,
                points="all",
                pointpos=0,
                marker=dict(size=2),
                jitter=0.6,
            ))
            violin_layout["xaxis"]["showticklabels"] = False
            violin_layout["xaxis_title"] = ""


        violin_layout["yaxis_title"] = feature
        fig.update_layout(violin_layout)

        return [fig]
    
        
    @staticmethod
    def get_callback_outputs():
        return [
            Output(component_id="violin-plot", component_property="figure"),
        ]
        
    # Inputs to Projection
    @staticmethod
    def get_callback_inputs():
        return {
            "feature": Input(component_id="violin-feature", component_property="value"),
            "groupby": Input(component_id="violin-groupby", component_property="value"),
        }

    @staticmethod
    def create_layout(dataset):
        var_names = [(f, f) for f in sorted(list(dataset.adata.var_names))]
        other = [(f, f.replace("_", " ").capitalize()) for f in sorted(list(dataset.adata.obs.columns)) if type(dataset.adata.obs[f][0]) != str and type(dataset.adata.obs[f][0]) != pd.CategoricalDtype]
        features = dict(other + var_names)


        groupbys = dict([(k, k.replace("_", " ").capitalize()) for k in dataset.adata.obs.columns if type(dataset.adata.obs[k][0]) == str or type(dataset.adata.obs[k][0]) == pd.CategoricalDtype])
        figure = html.Div(children=[
            html.Div(children=[ 
                # Features
                html.Div(children=[
                    html.Label("Feature"),
                    dcc.Dropdown(features, value=list(features.keys())[0], id="violin-feature", clearable=False),
                ], style={"flex": "1"}),
                # Groupby select
                html.Div(children=[
                    html.Label("Group By"),
                    dcc.Dropdown(groupbys, value=None, id="violin-groupby", clearable=True),
                ], style={"flex": "1"}),
            ], id="violin-select"),
            html.Div([
                dcc.Loading(
                    id="violin-projection", type="circle",
                    children=[html.Div(dcc.Graph(id="violin-plot"))],
                )
            ], id="violin-figure")
        ], id="violin")

        return figure