import dash
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

from cellbro.plots.Heatmap import Heatmap, heatmap_params, AddGenesFromList
import cellbro.plots.Projection as Projection
from cellbro.plots.Trimap import Trimap
from cellbro.plots.TSNE import TSNE
from cellbro.plots.UMAP import UMAP, SCVI_UMAP
from cellbro.plots.Violin import Violin
from cellbro.util.DashAction import DashAction
from cellbro.util.DashPage import DashPage
import cellbro.util.Components as Components

import scout

class PlotProjection(DashAction):
    def plot(self, params):
        color = params.pop("projection_color")
        obsm_layer = params.pop("projection_type")

        fig = scout.ply.projection(
            self.dataset.adata, obsm_layer=obsm_layer, hue=color,
            layout=Projection.projection_layout
        )

        return [fig]
        

    def setup_callbacks(self, app):
        outputs = [
            Output(component_id="projection-plot", component_property="figure"),
            # Output(component_id="projection-umap", component_property="style"),
            # Output(component_id="projection-tsne", component_property="style"),
            # Output(component_id="projection-trimap", component_property="style"),
            # Output(component_id="projection-scvi_umap", component_property="style"),
            # Output(component_id="projection-pca", component_property="style"),
        ]

        inputs = {
            # "submit": Input(
            #     component_id="projection-submit", component_property="n_clicks"
            # ),
            "projection_color": Input(
                component_id="projection-color", component_property="value"
            ),
            "projection_type": Input(
                component_id="projection-type", component_property="value"
            ),
        }

        states = {}
        for key in UMAP._params.keys():
            states[f"umap_{key}"] = State(
                component_id=f"projection-umap-{key}", component_property="value"
            )

        for key in SCVI_UMAP._params.keys():
            states[f"scvi_umap_{key}"] = State(
                component_id=f"projection-scvi_umap-{key}", component_property="value"
            )

        for key in TSNE._params.keys():
            states[f"tsne_{key}"] = State(
                component_id=f"projection-tsne-{key}", component_property="value"
            )

        for key in Trimap._params.keys():
            states[f"trimap_{key}"] = State(
                component_id=f"projection-trimap-{key}", component_property="value"
            )

        # for key in pca_params.keys():
        #     states[f"pca_{key}"] = State(
        #         component_id=f"projection-pca-{key}", component_property="value"
        #     )
        # Projection
        @app.dash_app.callback(output=outputs, inputs=inputs, state=states)
        def _(**kwargs):
            return self.plot(params=kwargs)

class PlotHeatmap(DashAction):
    def apply(self, params):
        return Heatmap(self.dataset, params).plot()

    def setup_callbacks(self, app):
        output = [
            Output(component_id="heatmap-plot", component_property="figure"),
            Output(component_id="heatmap-plot", component_property="style"),
        ]

        # Inputs to Projection
        inputs = {
            "submit": Input(
                component_id="heatmap-submit", component_property="n_clicks"
            ),
        }

        state = dict(
            [
                (
                    key,
                    State(component_id=f"heatmap-{key}", component_property="value"),
                )
                for key in heatmap_params.keys()
            ]
        )
        state["selected_genes"] = State("heatmap-selected-genes", "value")

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, **kwargs):
            return self.apply(params=kwargs)

class PlotViolin(DashAction):
    def apply(self, params):
        return Violin(self.dataset).plot(**params)

    def setup_callbacks(self, app):
        output = [
            Output(component_id="violin-plot", component_property="figure"),
        ]

        inputs = {
            "feature": Input(component_id="violin-feature", component_property="value"),
            "groupby": Input(component_id="violin-groupby", component_property="value"),
        }
        callbacks = dict(output=output, inputs=inputs)

        @app.dash_app.callback(**callbacks)
        def _(**kwargs):
            return self.apply(params=kwargs)


class UpdateAvailableProjectionTypes(DashAction):
    def setup_callbacks(self, app):
        @app.dash_app.callback(
            output=Output("projection-type", "options"),
            inputs={"_": Input("projection-type", "value")}
        )
        def _(_):
            # neighbors = self.dataset.get_neighbors()
            # available_projections = ["UMAP", "Trimap", "t-SNE", "PCA"]
            # if "neighbors_scvi" in neighbors:
            #     available_projections.append("SCVI-UMAP")
            available_projections = list(self.dataset.adata.obsm.keys())

            return available_projections

class CellsPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.cells", "Cells", "/cells", order)
        self.dataset = dataset
        self.actions.update(
            projection_action=PlotProjection(self.dataset),
            heatmap_action=PlotHeatmap(self.dataset),
            violin_action=PlotViolin(self.dataset),
            update_projection_types=UpdateAvailableProjectionTypes(self.dataset),
            heatmap_add_genes=AddGenesFromList(self.dataset),
        )

    def _params_layout(self):
        return [
            html.Div(
                children=Projection.Projection.get_layout(
                    UMAP),
                style={"display": "none"},
                id=f"projection-{UMAP.get_type().value}",
                className="param-class"
            ),
            html.Div(
                children=Projection.Projection.get_layout(
                    TSNE),
                style={"display": "none"},
                id=f"projection-{TSNE.get_type().value}",
                className="param-class"
            ),
            html.Div(
                children=Projection.Projection.get_layout(
                    Trimap),
                style={"display": "none"},
                id=f"projection-{Trimap.get_type().value}",
                className="param-class"
            ),
            html.Div(
                children=Projection.Projection.get_layout(
                    SCVI_UMAP),
                style={"display": "none"},
                id=f"projection-{SCVI_UMAP.get_type().value}",
                className="param-class"
            ),
        ]

    def create_layout(self) -> list:
        top_sidebar = Components.create_sidebar(
            id="cells-top-sidebar", class_name="top-sidebar",
            title="Projection Settings",
            params_children=self._params_layout(),
            btn_id="projection-submit", btn_text="Plot"
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
                                    id="projection-type", clearable=False,
                                ),
                            ],
                            className="param-column",
                        ),
                        # Projection Hue celect
                        html.Div(
                            children=[
                                html.Label("Color"),
                                dcc.Dropdown(
                                    self.dataset.adata.obs_keys()
                                    + list(self.dataset.adata.var_names),
                                    value=self.dataset.adata.obs_keys()[0],
                                    id="projection-color",
                                    clearable=False,
                                ),
                            ],
                            className="param-column",
                        ),
                    ],
                    id="projection-select",
                    className="main-select top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            id="loading-projection",
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id="projection-plot", className="main-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    id="projection-figure",
                    className="main-figure",
                ),
            ],
            className="main",
        )

        bottom_left_sidebar, bottom_figure = Heatmap.create_layout(self.dataset)
        violin_layout = Violin.create_layout(self.dataset)

        return [
            html.Div(
                id="top",
                className="top",
                children=[top_sidebar, main_figure, violin_layout],
            ),
            # html.Div(id="resizer", className="resizer"),
            html.Div(
                id="bottom",
                className="bottom",
                children=[bottom_left_sidebar, bottom_figure],
            ),
        ]
