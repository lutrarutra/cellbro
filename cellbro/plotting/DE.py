import plotly.graph_objects as go
import plotly.express as px
from dash import html, dcc, Input, Output, State
import dash_bootstrap_components as dbc
from plotly.subplots import make_subplots

from cellbro.util.Param import *

import scanpy as sc
import pandas as pd
import numpy as np

import scout

figure_layout = go.Layout(
    paper_bgcolor='white',
    plot_bgcolor='white',
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=0, b=0, l=0, r=0),
)

de_params = ParamsDict([
    Param(
        "method", "Method", default="t-test", type=str, description="",
        allowed_values=["logreg", "t-test", "wilcoxon", "t-test_overestim_var"]
    ),
    Param(
        "corr_method", "Correction Method", default="benjamini-hochberg", type=str, description="",
        allowed_values=["benjamini-hochberg", "t-test", "t-test_overestim_var", "wilcoxon"]
    ),
])

class DE():
    def __init__(self, dataset, params):
        self.dataset = dataset
        self.params = params

    def apply(self):
        print(self.params["submit"])
        if self.params["submit"] is None:
            groupby_options = DE._get_groupby_options(self.dataset)
            if len(groupby_options) == 0:
                return [[],None, [], None]
            
            groupby_selected = groupby_options[0]
            refs_options = DE._get_reference_options(self.dataset, groupby_selected)
            refs_selected = next(iter(self.dataset.adata.uns[f"rank_genes_{groupby_selected}"].keys()))
            return [groupby_options, groupby_selected, refs_options, refs_selected]

        key = f"rank_genes_{self.params['groupby']}"

        if key not in self.dataset.adata.uns.keys():
            scout.tl.rank_marker_genes(self.dataset.adata, groupby=self.params["groupby"])

        groupby_options = DE._get_groupby_options(self.dataset)
        refs_options = DE._get_reference_options(self.dataset, self.params["groupby"])
        
        groupby_selected = self.params["groupby"]
        refs_selected = next(iter(self.dataset.adata.uns[key].keys()))
        return [groupby_options, groupby_selected, refs_options, refs_selected]

    def de_volcano(self):
        key = f"rank_genes_{self.params['groupby']}"
        if self.params["reference"] is None:
            self.params["reference"] = list(self.dataset.adata.uns[key].keys())[0]
        fig = scout.ply.marker_volcano(
            self.dataset.adata.uns[key][self.params["reference"]],
            layout=figure_layout
        )
        return fig

    def plot(self):
        volcano = self.de_volcano()
        return [volcano]

    @staticmethod
    def get_apply_callbacks():
        outputs = [
            Output(component_id="volcano-groupby", component_property="options"),
            Output(component_id="volcano-groupby", component_property="value"),
            Output(component_id="volcano-reference", component_property="options"),
            Output(component_id="volcano-reference", component_property="value"),
        ]
        inputs = {
            "submit": Input(component_id="de-submit", component_property="n_clicks"),
        }
        states = {
            "groupby": State(component_id="de-groupby", component_property="value"),
        }
        for param in de_params.values():
            states[param.key] = State(component_id=f"de-{param.key}", component_property="value")

        return dict(output=outputs, inputs=inputs, state=states)

    @staticmethod
    def get_plot_callbacks():
        outputs = [
            Output(component_id="de-volcano-plot", component_property="figure"),
        ]
        inputs = {
            "groupby": Input(component_id="volcano-groupby", component_property="value"),
            "reference": Input(component_id="volcano-reference", component_property="value"),
        }
        states = {}
        return dict(output=outputs, inputs=inputs, state=states)

    @staticmethod
    def params_layout(dataset):
        cats = dataset.get_categoricals()
        divs = []

        divs.append(
            html.Div(children=[
                html.Label(
                    "Group By", className="param-label",
                ),
                html.Div([
                    dcc.Dropdown(
                        options=cats, value=cats[0], id=f"de-groupby", clearable=False
                    )
                ], className="param-select")
            ], className="param-row-stacked")
        )

        for key, param in de_params.items():
            divs.append(
                html.Div(children=[
                    html.Label(
                        param.name, className="param-label",
                    ),
                    html.Div([
                        dcc.Dropdown(
                        id=f"de-{key}", value=param.default, options=param.allowed_values, clearable=False
                        )
                    ], className="param-select"),
            ], className="param-row-stacked"))

        layout = html.Div(children=divs)
        return layout

    @staticmethod
    def _get_groupby_options(dataset):
        return [x[11:] for x in dataset.adata.uns.keys() if x.startswith("rank_genes_")]

    @staticmethod
    def _get_reference_options(dataset, groupby):
        if groupby is None:
            return []
        return dict([(x, f"{x} vs. Rest") for x in dataset.adata.uns[f"rank_genes_{groupby}"].keys()])

    @staticmethod
    def create_layout(dataset):
        top_sidebar = html.Div(children=[
            html.Div([
                html.H3("Differential Expression Settings"),
            ], className="sidebar-header"),
            dcc.Loading(type="circle", children=[
                html.Div(children=[
                    DE.params_layout(dataset),
                ], className="sidebar-parameters"),
                html.Div([
                    dbc.Button("Apply", color="primary", className="mr-1", id="de-submit"),
                ], className="sidebar-footer")
            ],),
        ], className="top-sidebar sidebar")

        groups = DE._get_groupby_options(dataset)
        refs = DE._get_reference_options(dataset, next(iter(groups), None))

        main = html.Div(children=[
            html.Div(children=[
                # Volcano Group By Select
                dcc.Loading(type="circle", children=[
                    html.Label("Volcano GroupBy"),
                    dcc.Dropdown(options=groups, value=next(iter(groups), None), id="volcano-groupby", clearable=False),
                ], className="param-column"),
                
                # Volcano Reference Select i.e. 'KO vs. Rest'
                dcc.Loading(type="circle", children=[
                    html.Label("Volcano Reference"),
                    dcc.Dropdown(options=refs, value=next(iter(refs), None), id="volcano-reference", clearable=False),
                ], className="param-column"),

            ], id="volcano-select", className="main-select top-parameters"),

            html.Div([
                dcc.Loading(
                    id="loading-de-volcano", type="circle",
                    children=[html.Div(dcc.Graph(id="de-volcano-plot", className="main-plot"))],
                )
            ], id="de-volcano-figure", className="main-figure")
        ], className="main")

        # bottom_sidebar = html.Div(children=[
        #     html.Div([
        #         html.H3("QC Violin Plots"),
        #     ], className="sidebar-header"),
        #     dcc.Loading(type="circle", children=[
        #         # html.Div(children=[
        #         # ], className="sidebar-parameters"),
        #         # html.Div([
        #         #     dbc.Button("Plot", color="primary", className="mr-1", id="qc-submit"),
        #         # ], className="sidebar-footer")
        #     ],)
        # ], id="qc-violin-sidebar", className="bottom-sidebar sidebar")


        # secondary_figure = html.Div(children=[
        #     html.Div([
        #         dcc.Loading(
        #             id="loading-dispersion", type="circle",
        #             children=[html.Div(dcc.Graph(id="dispersion-plot", className="secondary-plot"))],
        #         )
        #     ], id="dispersion-figure", className="secondary-figure")
        # ], className="secondary")

        # bottom_figure = html.Div(children=[
        #     dcc.Loading(
        #         id="loading-qc-violin", className="loading-bottom", type="circle",
        #         children=[html.Div(dcc.Graph(id="qc-violin-plot", className="bottom-plot"))],
        #     )
        # ], id="qc-violin-figure", className="bottom-figure")

        return top_sidebar, main