import dash
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, State, dcc, html

from cellbro.plots.Heatmap import Heatmap, heatmap_params
from cellbro.plots.Projection import Projection, ProjectionType
from cellbro.plots.Trimap import Trimap
from cellbro.plots.TSNE import TSNE
from cellbro.plots.UMAP import UMAP
from cellbro.plots.Violin import Violin
from cellbro.util.DashAction import DashAction
from cellbro.util.DashPage import DashPage


def parse_params(params):
    projection_type = params.pop("projection_type")
    projection_color = params.pop("projection_color")
    key = None
    if projection_type == "UMAP":
        key = "umap"
    elif projection_type == "t-SNE":
        key = "tsne"
    else:
        key = "trimap"

    projection_params = dict(
        [
            (param_key.replace(f"{key}_", ""), params[param_key])
            for param_key in params.keys()
            if param_key.startswith(f"{key}_")
        ]
    )
    return projection_type, dict(color=projection_color, params=projection_params)


def projection_factory(dataset, params) -> Projection:
    projection_type, kwargs = parse_params(params)
    if projection_type == "UMAP":
        return UMAP(dataset=dataset, **kwargs)
    elif projection_type == "t-SNE":
        return TSNE(dataset=dataset, **kwargs)
    # elif projection_type == "Trimap":
    else:
        return Trimap(dataset=dataset, **kwargs)
    # elif projection_type == "PCA":
    #     return PCA(dataset=dataset, color=color, params=params)


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
                # {"display": "none"},
            ]

        elif projection.get_type() == ProjectionType.TSNE:
            return [
                projection.plot(),
                {"display": "none"},
                {"display": "block"},
                {"display": "none"},
                # {"display": "none"},
            ]

        elif projection.get_type() == ProjectionType.TRIMAP:
            return [
                projection.plot(),
                {"display": "none"},
                {"display": "none"},
                {"display": "block"},
                # {"display": "none"},
            ]
        assert False, "Invalid projection type"

    def setup_callbacks(self, dash_app):
        outputs = [
            Output(component_id="projection-plot", component_property="figure"),
            Output(component_id="projection-umap", component_property="style"),
            Output(component_id="projection-tsne", component_property="style"),
            Output(component_id="projection-trimap", component_property="style"),
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

        for key in TSNE._params.keys():
            states[f"tsne_{key}"] = State(
                component_id=f"projection-tsne-{key}", component_property="value"
            )

        for key in Trimap._params.keys():
            states[f"trimap_{key}"] = State(
                component_id=f"projection-trimap-{key}", component_property="value"
            )
        callbacks = dict(output=outputs, inputs=inputs, state=states)

        # for key in pca_params.keys():
        #     states[f"pca_{key}"] = State(
        #         component_id=f"projection-pca-{key}", component_property="value"
        #     )
        # Projection
        @dash_app.callback(**callbacks)
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


class CellsPage(DashPage):
    def __init__(self, dataset, dash_app):
        super().__init__("pages.cells", "Cells", "/cells", 2)
        self.dataset = dataset
        self.actions = dict(
            projection_action=PlotProjection(self.dataset),
            heatmap_action=PlotHeatmap(self.dataset),
            violin_action=PlotViolin(self.dataset),
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
                                Projection.get_layout(UMAP),
                                Projection.get_layout(TSNE),
                                Projection.get_layout(Trimap),
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

        main_figure = html.Div(
            children=[
                html.Div(
                    children=[
                        # Projection type celect
                        html.Div(
                            children=[
                                html.Label("Projection Type"),
                                dcc.Dropdown(
                                    ["UMAP", "Trimap", "t-SNE", "PCA"],
                                    value="UMAP",
                                    id="projection-type",
                                    clearable=False,
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
