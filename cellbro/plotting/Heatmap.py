import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import Input, Output, State, dcc, html

from cellbro.util.Param import *

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
            default="log1p",
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
            default="RdBu_r",
            type=list,
            description="",
            allowed_values={
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
        if self.params["layer"] == "log1p":
            z = self.dataset.adata[:, self.selected_genes].X.toarray()
        else:
            z = (
                self.dataset.adata[:, self.selected_genes]
                .layers[self.params["layer"]]
                .toarray()
            )
        fig = px.imshow(
            z.T,
            y=self.selected_genes,
            aspect="auto",
            color_continuous_scale=self.params["colormap"],
            color_continuous_midpoint=0 if "centered" in self.params["layer"] else None,
        )
        fig.update_layout(heatmap_layout)
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
            "height": f"{z.shape[1] * 20 + 50}px",
        }
        return fig, style

    @staticmethod
    def get_callback_outputs():
        return [
            Output(component_id="heatmap-plot", component_property="figure"),
            Output(component_id="heatmap-plot", component_property="style"),
        ]

    # Inputs to Projection
    @staticmethod
    def get_callback_inputs():
        return {
            "submit": Input(
                component_id="heatmap-submit", component_property="n_clicks"
            ),
        }

    @staticmethod
    def get_callback_states():
        states = dict(
            [
                (key, State(component_id=f"heatmap-{key}", component_property="value"))
                for key in heatmap_params.keys()
            ]
        )
        return states

    @staticmethod
    def params_layout():
        divs = []
        for key, param in heatmap_params.items():
            if param.type == list:
                divs.append(
                    html.Div(
                        children=[
                            html.Label(
                                param.name,
                                className="param-label",
                            ),
                            html.Div(
                                [
                                    dcc.Dropdown(
                                        options=param.allowed_values,
                                        value=param.default,
                                        id=f"heatmap-{key}",
                                        clearable=False,
                                    )
                                ],
                                className="param-select",
                            ),
                        ],
                        className="param-row-stacked",
                    )
                )
            else:
                divs.append(
                    html.Div(
                        [
                            html.Label(
                                param.name,
                                className="param-label",
                            ),
                            dcc.Input(
                                id=f"heatmap-{key}",
                                type=param.input_type,
                                value=param.value,
                                step=param.step if param.step != None else 0.1,
                                className="param-input",
                            ),
                        ],
                        className="param-row",
                    )
                )
        return html.Div(children=divs)

    @staticmethod
    def create_layout(dataset):
        sidebar = html.Div(
            children=[
                html.Div(
                    [
                        html.H3("Heatmap Settings"),
                    ],
                    id="heatmap-header",
                    className="sidebar-header",
                ),
                dcc.Loading(
                    type="circle",
                    children=[
                        html.Div(
                            children=[
                                Heatmap.params_layout(),
                            ],
                            className="sidebar-parameters",
                        ),
                        html.Div(
                            [
                                dbc.Button(
                                    "Plot",
                                    color="primary",
                                    className="mr-1",
                                    id="heatmap-submit",
                                ),
                            ],
                            id="heatmap-footer",
                            className="sidebar-footer",
                        ),
                    ],
                ),
            ],
            id="heatmap-sidebar",
            className="bottom-sidebar sidebar",
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
