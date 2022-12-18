import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

from cellbro.util.Param import *
import cellbro.util.Components as Components
from cellbro.util.DashAction import DashAction

import scout

heatmap_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=10, b=10, l=10, r=10),
)

heatmap_params = ParamsDict(
    [
        Param(
            key="layer",
            name="Layer",
            default="logcentered",
            type=list,
            description="",
            allowed_values={
                "log1p": "log1p",
                "counts": "Counts",
                "ncounts": "Normalized Counts",
                "centered": "Centered Counts",
                "logcentered": "Log Centered",
            },
        ),
        Param(
            key="colormap",
            name="Colormap",
            default="seismic",
            type=list,
            description="",
            allowed_values={
                "seismic": "Seismic (for centered)",
                "RdBu_r": "B-W-R",
                "viridis": "Viridis",
                "plasma": "Plasma",
                "inferno": "Inferno",
                "magma": "Magma",
                "cividis": "Cividis",
            },
        ),
    ]
)

class AddGenesFromList(DashAction):
    def setup_callbacks(self, app):
        output = Output("heatmap-selected-genes", "value")
        inputs = [Input("heatmap-selected-genelists", "value")]
        state = [State("heatmap-selected-genes", "value")]
        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(gene_lists, selected_genes):
            if gene_lists is None:
                raise PreventUpdate

            if selected_genes is None:
                selected_genes = []

            for gene_list in gene_lists:
                selected_genes.extend(self.dataset.get_genes_from_list(gene_list))

            return list(set(selected_genes))

class Heatmap:
    def __init__(self, dataset, params):
        self.dataset = dataset
        self.params = params
        selected_genes = params.pop("selected_genes")
        if selected_genes and len(selected_genes) > 0:
            self.selected_genes = selected_genes
        else:
            self.selected_genes = self.dataset.adata.var_names[:50].tolist()

    def plot(self):
        fig = scout.ply.heatmap(
            adata=self.dataset.adata, var_names=self.selected_genes,
            categoricals=["leiden"], layer=self.params["layer"], cluster_cells=False, layout=dict()
        )

        style = {
            # "width": f"{int(z.shape[0]/2)+100}px",
            "height": fig.layout["height"],
        }
        return fig, style

    @staticmethod
    def create_layout(dataset):
        sidebar = Components.create_sidebar(
            id="cells-bot-sidebar", class_name="bot-sidebar",
            title="Heatmap Settings",
            params_children=Heatmap._params_layout(dataset),
            btn_id="heatmap-submit", btn_text="Plot"
        )

        figure = html.Div(
            children=[
                html.Div(
                    [
                        dcc.Loading(
                            id="loading-heatmap",
                            className="loading-bottom",
                            type="circle",
                            children=[
                                html.Div(
                                    [
                                        dcc.Graph(
                                            id="heatmap-plot", className="bottom-plot"
                                        )
                                    ],
                                    id="heatmap-plot-colorbar",
                                )
                            ],
                        )
                    ],
                    id="heatmap-figure",
                    className="bottom-figure",
                )
            ]
        )

        return sidebar, figure

    @staticmethod
    def _params_layout(dataset):
        genes = sorted(dataset.adata.var_names.tolist())
        gene_lists = sorted(dataset.get_gene_lists())

        divs = [
            html.Div([
                html.Label(
                    "Show Genes",
                    className="param-label",
                ),
                html.Div([
                    dcc.Dropdown(
                        options=genes, value=None, id="heatmap-selected-genes", clearable=True,
                        placeholder="Select Genes", multi=True,
                        style={"width": "100%"}
                    ),
                    dcc.Dropdown(
                        options=gene_lists, value=None, id="heatmap-selected-genelists", clearable=True,
                        placeholder="Select Gene Lists", multi=True,
                        style={"width": "100%"}
                    ),
                ], style={"display":"flex", "gap":"10px"})
            ], className="param-row-stacked")
        ]
        divs.extend(Components.params_layout(heatmap_params, "heatmap"))

        return divs

