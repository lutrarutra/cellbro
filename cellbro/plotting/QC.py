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
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=0, b=0, l=0, r=0),
)

qc_params = ParamsDict([
    Param(key="pct_counts_mt", name="Max MT%", default=5.0, type=float, description="", step=1),
    Param(key="min_genes", name="Min. Genes (per cell)", default=200, type=int, description="", step=100),
    Param(key="min_cells", name="Min. Cells (per gene)", default=10, type=int, description="", step=10),
])

class QC():
    def __init__(self, dataset, params):
        self.dataset = dataset
        self.params = params

    def filter(self, submit):
        # Makes sure that filtering is not done on initial load
        if submit is None:
            return list(self.params.values())
        print("Filtering...")
        self.dataset.adata = self.dataset.adata[self.dataset.adata.obs.pct_counts_mt < self.params["pct_counts_mt"], :].copy()
        sc.pp.filter_cells(self.dataset.adata, min_genes=self.params["min_genes"])
        sc.pp.filter_genes(self.dataset.adata, min_cells=self.params["min_cells"])
        return list(self.params.values())

    def violin_plot(self):
        violins = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
        fig = make_subplots(rows=1, cols=3)
        x = self.dataset.adata.obs[violins]
        df = pd.DataFrame(
            x, columns=violins
        )

        for i, feature in enumerate(violins):
            fig.add_trace(
                    go.Violin(
                        name=feature.replace("_", " ").title(),
                        y=df[feature],
                        box_visible=True,
                        points="all",
                        pointpos=0,
                        marker=dict(size=2),
                        jitter=0.6,
                ), row=1, col=i+1
            )

        fig.update_layout(figure_layout)
        fig.update_layout(showlegend=False)

        return fig
    
    def mt_plot(self):
        color = (self.dataset.adata.obs["pct_counts_mt"] > self.params["pct_counts_mt"]).values
        
        cmap = ["#d3d3d3", sc.pl.palettes.default_20[0]] if color[0] else [sc.pl.palettes.default_20[0], "#d3d3d3"]

        fig = px.scatter(
            self.dataset.adata.obs,
            x="total_counts",
            y="pct_counts_mt",
            color=color,
            color_discrete_sequence=cmap,
        )
        fig.update_layout(figure_layout)
        fig.update_layout(xaxis_title="Total Counts", yaxis_title="MT%", legend_title_text=f"MT > {self.params['pct_counts_mt']:.1f}%")
        fig.add_hline(y=self.params["pct_counts_mt"], line_width=1, line_dash="dash", line_color=sc.pl.palettes.default_20[3])
        return fig

    def plot(self):
        main_figure = self.mt_plot()
        bottom_figure = self.violin_plot()

        return [main_figure, bottom_figure]

    @staticmethod
    def get_filtering_callbacks():
        outputs = []
        for param in qc_params.values():
            outputs.append(Output(component_id=f"qc-{param.key}", component_property="value"))

        inputs = {
            "submit": Input(component_id="filtering-submit", component_property="n_clicks"),
        }
        states = {}
        for param in qc_params.values():
            states[param.key] = State(component_id=f"qc-{param.key}", component_property="value")
        return outputs, inputs, states
    
    @staticmethod
    def get_callback_outputs():
        return [
            Output(component_id="mt-plot", component_property="figure"),
            Output(component_id="qc-violin-plot", component_property="figure"),
        ]
        
    # Inputs to Projection
    @staticmethod
    def get_callback_inputs():
        inputs = {}
        for param in qc_params.values():
            inputs[param.key] = Input(component_id=f"qc-{param.key}", component_property="value")
        return inputs

    @staticmethod
    def create_layout(dataset):
        top_sidebar = html.Div(children=[
            html.Div([
                html.H3("Filtering Settings"),
            ], className="top-header"),
            dcc.Loading(type="circle", children=[
                html.Div(children=[
                    QC.params_layout(),
                ], className="top-parameters"),
                html.Div([
                    dbc.Button("Filter", color="primary", className="mr-1", id="filtering-submit"),
                ], className="top-footer")
            ],),
        ], className="top-sidebar")


        bottom_sidebar = html.Div(children=[
            dcc.Loading(type="circle", children=[
                html.Div([
                    dbc.Button("Plot", color="primary", className="mr-1", id="qc-submit"),
                ], id="bottom-sidebar-footer")
            ],)
        ], id="qc-violin-sidebar", className="bottom-sidebar")

        main_figure = html.Div(children=[
            html.Div([
                dcc.Loading(
                    id="loading-mt", type="circle",
                    children=[html.Div(dcc.Graph(id="mt-plot", className="main-plot"))],
                )
            ], id="mt-figure", className="main-figure")
        ], className="main")

        bottom_figure = html.Div(children=[
            dcc.Loading(
                id="loading-qc-violin", className="loading-bottom", type="circle",
                children=[html.Div(dcc.Graph(id="qc-violin-plot", className="bottom-plot"))],
            )
        ], id="qc-violin-figure", className="bottom-figure")

        return top_sidebar, main_figure, bottom_sidebar, bottom_figure

    @staticmethod
    def params_layout():
        layout = html.Div(children=[
            html.Div(children=[
                html.Label(
                    param.name, className="param-label",
                ),
                dcc.Input(
                    id=f"qc-{key}", type=param.input_type, value=param.value,
                    step=param.step if param.step != None else 0.1, className="param-input",
                ),
            ], className="param-row") for key, param in qc_params.items()
        ])
        return layout