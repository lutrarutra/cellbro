import dash
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, State, dcc, html

from cellbro.plots.Heatmap import Heatmap, heatmap_params
from cellbro.plots.Projection import Projection, ProjectionType, parse_params
from cellbro.plots.Trimap import Trimap
from cellbro.plots.TSNE import TSNE
from cellbro.plots.UMAP import UMAP, SCVI_UMAP
from cellbro.plots.Violin import Violin
from cellbro.util.DashAction import DashAction
from cellbro.util.DashPage import DashPage


def projection_factory(dataset, params) -> Projection:
    projection_type, kwargs = parse_params(params)
    if projection_type == "UMAP":
        return UMAP(dataset=dataset, **kwargs)
    elif projection_type == "t-SNE":
        return TSNE(dataset=dataset, **kwargs)
    # elif projection_type == "Trimap":
    elif projection_type == "SCVI-UMAP":
        return SCVI_UMAP(dataset=dataset, **kwargs)
    elif projection_type == "Trimap":
        return Trimap(dataset=dataset, **kwargs)
    else:
        assert False, "Invalid projection type"

class PlotProjection(DashAction):
    def apply(self, params):
        projection = projection_factory(self.dataset, params)
        projection.apply()
        if projection.get_type() == ProjectionType.UMAP:
            return [
                projection.plot(),
                {"display": "block"},
                {"display": "none"},
                {"display": "none"},
                {"display": "none"},
            ]

        elif projection.get_type() == ProjectionType.TSNE:
            return [
                projection.plot(),
                {"display": "none"},
                {"display": "block"},
                {"display": "none"},
                {"display": "none"},
            ]

        elif projection.get_type() == ProjectionType.TRIMAP:
            return [
                projection.plot(),
                {"display": "none"},
                {"display": "none"},
                {"display": "block"},
                {"display": "none"},
            ]
        elif projection.get_type() == ProjectionType.SCVI_UMAP:
            return [
                projection.plot(),
                {"display": "none"},
                {"display": "none"},
                {"display": "none"},
                {"display": "block"},
            ]
        assert False, "Invalid projection type"

    def setup_callbacks(self, dash_app):
        outputs = [
            Output(component_id="projection-plot", component_property="figure"),
            Output(component_id="projection-umap", component_property="style"),
            Output(component_id="projection-tsne", component_property="style"),
            Output(component_id="projection-trimap", component_property="style"),
            Output(component_id="projection-scvi_umap", component_property="style"),
            # Output(component_id="projection-pca", component_property="style"),
        ]

        inputs = {
            "submit": Input(
                component_id="projection-submit", component_property="n_clicks"
            ),
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
        @dash_app.callback(output=outputs, inputs=inputs, state=states)
        def _(**kwargs):
            return self.apply(params=kwargs)


class PlotHeatmap(DashAction):
    def apply(self, params):
        return Heatmap(self.dataset, params).plot()

    def setup_callbacks(self, dash_app):
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
        callbacks = dict(output=output, inputs=inputs, state=state)

        @dash_app.callback(**callbacks)
        def _(submit, **kwargs):
            return self.apply(params=kwargs)

class PlotViolin(DashAction):
    def apply(self, params):
        return Violin(self.dataset).plot(**params)

    def setup_callbacks(self, dash_app):
        output = [
            Output(component_id="violin-plot", component_property="figure"),
        ]

        inputs = {
            "feature": Input(component_id="violin-feature", component_property="value"),
            "groupby": Input(component_id="violin-groupby", component_property="value"),
        }
        callbacks = dict(output=output, inputs=inputs)

        @dash_app.callback(**callbacks)
        def _(**kwargs):
            return self.apply(params=kwargs)


class UpdateAvailableProjectionTypes(DashAction):
    def setup_callbacks(self, dash_app):
        @dash_app.callback(
            output=Output("projection-type", "options"),
            inputs={"_": Input("projection-type", "value")}
        )
        def _(_):
            neighbors = self.dataset.get_neighbors()
            available_projections = ["UMAP", "Trimap", "t-SNE", "PCA"]
            if "neighbors_scvi" in neighbors:
                available_projections.append("SCVI-UMAP")

            return available_projections

class CellsPage(DashPage):
    def __init__(self, dataset, dash_app, order):
        super().__init__("pages.cells", "Cells", "/cells", order)
        self.dataset = dataset
        self.actions = dict(
            projection_action=PlotProjection(self.dataset),
            heatmap_action=PlotHeatmap(self.dataset),
            violin_action=PlotViolin(self.dataset),
            update_projection_types=UpdateAvailableProjectionTypes(self.dataset),
        )
        self.setup_callbacks(dash_app)

    def create_layout(self) -> list:
        left_sidebar = html.Div(
            children=[
                html.Div(
                    [
                        html.H3("Projection Settings"),
                    ],
                    className="sidebar-header",
                ),
                dcc.Loading(
                    type="circle",
                    children=[
                        html.Div(
                            children=[
                                html.Div(
                                    children=Projection.get_layout(UMAP),
                                    style={"display": "none"},
                                    id=f"projection-{UMAP.get_type().value}",
                                    className="param-class"
                                ),
                                html.Div(
                                    children=Projection.get_layout(TSNE),
                                    style={"display": "none"},
                                    id=f"projection-{TSNE.get_type().value}",
                                    className="param-class"
                                ),
                                html.Div(
                                    children=Projection.get_layout(Trimap),
                                    style={"display": "none"},
                                    id=f"projection-{Trimap.get_type().value}",
                                    className="param-class"
                                ),
                                html.Div(
                                    children=Projection.get_layout(SCVI_UMAP),
                                    style={"display": "none"},
                                    id=f"projection-{SCVI_UMAP.get_type().value}",
                                    className="param-class"
                                ),
                            ],
                            className="sidebar-parameters",
                        ),
                        html.Div(
                            [
                                dbc.Button(
                                    "Plot",
                                    color="primary",
                                    className="mr-1",
                                    id="projection-submit",
                                ),
                            ],
                            className="sidebar-footer",
                        ),
                    ],
                ),
            ],
            className="top-sidebar sidebar",
        )
        neighbors = self.dataset.get_neighbors()
        available_projections = ["UMAP", "Trimap", "t-SNE", "PCA"]
        if "neighbors_scvi" in neighbors:
            available_projections.append("SCVI-UMAP")

        main_figure = html.Div(
            children=[
                html.Div(
                    children=[
                        # Projection type celect
                        html.Div(
                            children=[
                                html.Label("Projection Type"),
                                dcc.Dropdown(
                                    available_projections, value="UMAP",
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
                                    + self.dataset.adata.var_names.tolist(),
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
                children=[left_sidebar, main_figure, violin_layout],
            ),
            # html.Div(id="resizer", className="resizer"),
            html.Div(
                id="bottom",
                className="bottom",
                children=[bottom_left_sidebar, bottom_figure],
            ),
        ]
