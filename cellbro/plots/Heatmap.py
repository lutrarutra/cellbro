import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import Input, Output, State, dcc, html

from cellbro.util.Param import *
import cellbro.util.Components as Components

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

class Heatmap:
    def __init__(self, dataset, params):
        self.dataset = dataset
        self.params = params
        self.selected_genes = dataset.adata.var_names[:50]

    def plot(self):
        # if self.params["layer"] == "log1p":
        #     z = self.dataset.adata[:, self.selected_genes].X
        #     if isinstance(z, scipy.sparse.csr_matrix):
        #         z = z.toarray()
        # else:
        #     z = (
        #         self.dataset.adata[:, self.selected_genes]
        #         .layers[self.params["layer"]]
        #     )
        #     if isinstance(z, scipy.sparse.csr_matrix):
        #         z = z.toarray()
                
        # fig = px.imshow(
        #     z.T,
        #     y=self.selected_genes,
        #     aspect="auto",
        #     color_continuous_scale=self.params["colormap"],
        #     color_continuous_midpoint=0 if "centered" in self.params["layer"] else None,
        # )
        # fig.update_layout(heatmap_layout)
        fig = scout.ply.heatmap(
            adata=self.dataset.adata, var_names=self.selected_genes,
            categoricals=["leiden"], layer=self.params["layer"], cluster_cells=False, layout=dict()
        )
        # mn, mx = z.min(), z.max()
        # colorbar = px.scatter(
        #     x=[0,0], y=[0,0], color=[mn, mx],
        #     color_continuous_scale=self.params["colormap"],
        #     color_continuous_midpoint=0 if "centered" in self.params["layer"] else None
        # )
        # colorbar.update_layout(heatmap_layout)
        # fig.update_coloraxes(showscale=False)

        # fig.update_layout(xaxis=dict(rangeslider=dict(visible=True, yaxis=dict(rangemode="fixed"), thickness=0.01), type="linear"))
        style = {
            # "width": f"{int(z.shape[0]/2)+100}px",
            "height": f"{len(self.selected_genes) * 20 + 50}px",
        }
        return fig, style

    @staticmethod
    def create_layout(dataset):
        sidebar = Components.create_sidebar(
            id="cells-bot-sidebar", class_name="bot-sidebar",
            title="Heatmap Settings",
            params_children=Components.params_layout(heatmap_params, "heatmap"),
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
