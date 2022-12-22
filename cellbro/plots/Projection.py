import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import Input, Output, State, dcc, html, ctx

from cellbro.plots.DashFigure import DashFigure
from cellbro.util.DashAction import DashAction
from cellbro.plots.UMAP import UMAP, SCVI_UMAP
from cellbro.plots.TSNE import TSNE
from cellbro.plots.Trimap import Trimap
import cellbro.util.Components as Components

import scout


projection_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
    margin=dict(t=10, b=10, l=10, r=10),
)

# class ProjectionType(Enum):
#     UMAP = "umap"
#     TRIMAP = "trimap"
#     TSNE = "tsne"
#     MDE = "mde"
#     SCVI_UMAP = "scvi_umap"
    # SCVI_MDE = "scvi_mde"
    # PCA = "PCA"

class PlotProjection(DashAction):
    def plot(self, color, obsm_layer):
        fig = scout.ply.projection(
            self.dataset.adata, obsm_layer=obsm_layer, hue=color,
            layout=projection_layout
        )
        return fig

    def apply(self, projection_type, params):
        if projection_type == "UMAP":
            projection = UMAP(self.dataset, UMAP.parse_params(params))
        elif projection_type == "t-SNE":
            projection = TSNE(self.dataset, TSNE.parse_params(params))
        elif projection_type == "Trimap":
            projection = Trimap(self.dataset, Trimap.parse_params(params))

        return projection.apply()

    def setup_callbacks(self, app):
        outputs = [
            Output(component_id=f"{self.page_id_prefix}-projection-plot", component_property="figure"),
            Output(f"{self.page_id_prefix}-projection-plot-type", "options"),
            Output(component_id=f"{self.page_id_prefix}-projection-plot-type", component_property="value"),
        ]

        inputs = {
            "projection_submit": Input(
                component_id=f"{self.page_id_prefix}-projection-submit", component_property="n_clicks"
            ),
            "color": Input(
                component_id=f"{self.page_id_prefix}-projection-color", component_property="value"
            ),
            "obsm_layer": Input(
                component_id=f"{self.page_id_prefix}-projection-plot-type", component_property="value"
            ),
        }

        states = {}
        states["projection_type"] = State(component_id=f"{self.page_id_prefix}-projection-type-select", component_property="value")
        for key in UMAP._params.keys():
            states[f"umap_{key}"] = State(
                component_id=f"{self.page_id_prefix}-projection-umap-{key}", component_property="value"
            )

        for key in SCVI_UMAP._params.keys():
            states[f"scvi_umap_{key}"] = State(
                component_id=f"{self.page_id_prefix}-projection-scvi_umap-{key}", component_property="value"
            )

        for key in TSNE._params.keys():
            states[f"tsne_{key}"] = State(
                component_id=f"{self.page_id_prefix}-projection-tsne-{key}", component_property="value"
            )

        for key in Trimap._params.keys():
            states[f"trimap_{key}"] = State(
                component_id=f"{self.page_id_prefix}-projection-trimap-{key}", component_property="value"
            )

        @app.dash_app.callback(output=outputs, inputs=inputs, state=states)
        def _(projection_submit, color, obsm_layer, projection_type, **kwargs):
            if ctx.triggered_id == f"{self.page_id_prefix}-projection-submit":
                if projection_submit is not None:
                    obsm_layer = self.apply(projection_type, params=kwargs)

            return (self.plot(color=color, obsm_layer=obsm_layer),
                    list(self.dataset.adata.obsm.keys()), obsm_layer)


class SelectProjectionType(DashAction):
    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-projection-umap", component_property="style"),
            Output(component_id=f"{self.page_id_prefix}-projection-tsne", component_property="style"),
            Output(component_id=f"{self.page_id_prefix}-projection-trimap", component_property="style"),
        ]
        inputs = [
            Input(component_id=f"{self.page_id_prefix}-projection-type-select", component_property="value"),
        ]

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(projection_type):
            if projection_type == "UMAP":
                return {"display": "block"}, {"display": "none"}, {"display": "none"}

            if projection_type == "t-SNE":
                return {"display": "none"}, {"display": "block"}, {"display": "none"}

            if projection_type == "Trimap":
                return {"display": "none"}, {"display": "none"}, {"display": "block"}

class Projection(DashFigure):
    def __init__(self, dataset, page_id_prefix):
        super().__init__(dataset, page_id_prefix)
        self.actions.update(
            projection_action=PlotProjection(self.dataset, self.page_id_prefix),
            select_projection_type=SelectProjectionType(self.dataset, self.page_id_prefix),
        )

    def _params_layout(self):
        return [
            html.Div([
                html.Div(
                    children=[
                        html.Label(
                            "Projection Type",
                            className="param-label",
                        ),
                        dcc.Dropdown(
                            id=f"{self.page_id_prefix}-projection-type-select",
                            options=["UMAP", "t-SNE", "Trimap"],
                            value="UMAP",
                            clearable=False,
                        )
                    ],
                    className="param-row-stacked",
                )
            ]),
            html.Div(
                children=UMAP.get_layout(self.page_id_prefix),
                style={"display": "none"},
                id=f"{self.page_id_prefix}-projection-{UMAP.get_key()}",
                className="param-class"
            ),
            html.Div(
                children=TSNE.get_layout(self.page_id_prefix),
                style={"display": "none"},
                id=f"{self.page_id_prefix}-projection-{TSNE.get_key()}",
                className="param-class"
            ),
            html.Div(
                children=Trimap.get_layout(self.page_id_prefix),
                style={"display": "none"},
                id=f"{self.page_id_prefix}-projection-{Trimap.get_key()}",
                className="param-class"
            ),
            html.Div(
                children=SCVI_UMAP.get_layout(self.page_id_prefix),
                style={"display": "none"},
                id=f"{self.page_id_prefix}-projection-{SCVI_UMAP.get_key()}",
                className="param-class"
            ),
        ]

    def create_layout(self) -> list:
        top_sidebar = Components.create_sidebar(
            id=f"{self.page_id_prefix}-top-sidebar", class_name="top-sidebar",
            title="Projection Settings",
            params_children=self._params_layout(),
            btn_id=f"{self.page_id_prefix}-projection-submit", btn_text="Apply Projection"
        )

        available_projections = list(self.dataset.adata.obsm.keys())

        main_figure = html.Div(
            children=[
                html.Div(
                    children=[
                        # Projection type celect
                        html.Div(
                            children=[
                                html.Label("Projection Type"),
                                dcc.Dropdown(
                                    available_projections, value=available_projections[0],
                                    id=f"{self.page_id_prefix}-projection-plot-type", clearable=False,
                                ),
                            ],
                            className="param-column",
                        ),
                        # Projection Hue celect
                        html.Div(
                            children=[
                                html.Label("Color"),
                                dcc.Dropdown(
                                    self.dataset.get_obs_features(),
                                    value=self.dataset.adata.obs_keys()[0],
                                    id=f"{self.page_id_prefix}-projection-color",
                                    clearable=False,
                                ),
                            ],
                            className="param-column",
                        ),
                    ],
                    id=f"{self.page_id_prefix}-projection-select",
                    className="main-select top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id=f"{self.page_id_prefix}-projection-plot", className="main-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    id=f"{self.page_id_prefix}-projection-figure",
                    className="main-figure",
                ),
            ],
            className="main",
        )

        return top_sidebar, main_figure


