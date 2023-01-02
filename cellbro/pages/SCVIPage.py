import dash
from dash import Input, Output, State, dcc, html, ctx
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from ..plots import SCVI
from ..components import components
from ..plots.projection import Projection
from ..components.DashPage import DashPage
from ..util.DashAction import DashAction
from ..plots.projection.UMAP import SCVI_UMAP

import scout

class FitAction(DashAction):
    def apply(self, params: dict):
        setup_params = dict(
            [(key[6:], value) for key, value in params.items() if key.startswith("setup_")]
        )
        model_params = dict(
            [(key[6:], value) for key, value in params.items() if key.startswith("model_")]
        )
        train_params = dict(
            [(key[6:], value) for key, value in params.items() if key.startswith("train_")]
        )
        scvi_plots.setup(self.dataset, setup_params)
        scvi_plots.fit(self.dataset, model_params, train_params)

    def apply_projection(self, params):
        projection = SCVI_UMAP(self.dataset, SCVI_UMAP.parse_params(params))
        return projection.apply()

    def plot(self, color, obsm_layer):        
        fig = scout.ply.projection(
            self.dataset.adata, obsm_layer=obsm_layer, hue=color,
            layout=Projection.projection_layout
        )

        return fig

    def setup_callbacks(self, app):
        output = [
            Output(f"{self.page_id}-projection-plot", "figure"),
            Output(f"{self.page_id}-projection-plot-type", "options"),
            Output(component_id=f"{self.page_id}-projection-plot-type", component_property="value"),
        ]
        inputs = {
            "submit": Input(f"{self.page_id}-fit-submit", "n_clicks"),
            "color": Input(
                component_id=f"{self.page_id}-projection-color", component_property="value"
            ),
            "obsm_layer": Input(
                component_id=f"{self.page_id}-projection-plot-type", component_property="value"
            ),
        }
        state = {
            "setup_continuous_covariate_keys": Input(f"{self.page_id}-continuous_covariate_keys", "value"),
            "setup_categorical_covariate_keys": Input(f"{self.page_id}-categorical_covariate_keys", "value"),
            "setup_batch_key": Input(f"{self.page_id}-batch_key", "value"),
        }
        for key, param in SCVI.scvi_model_params.items():
            state[f"model_{param.key}"] = State(f"{self.page_id}-model-{param.key}", "value")

        for key, param in SCVI.scvi_train_params.items():
            state[f"train_{param.key}"] = State(f"{self.page_id}-train-{param.key}", "value")

        for key in SCVI_UMAP._params.keys():
            state[f"scvi_umap_{key}"] = State(
                component_id=f"{self.page_id}-projection-scvi_umap-{key}", component_property="value"
            )

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, color, obsm_layer, **kwargs):
            if ctx.triggered_id == f"{self.page_id}-fit-submit":
                self.apply(params=kwargs)
                obsm_layer = self.apply_projection(params=kwargs)

            layers = self.dataset.get_scvi_projections()
            if len(layers) == 0:
                raise PreventUpdate

            if obsm_layer is None:
                obsm_layer = layers[0]

            fig = self.plot(color, obsm_layer)
            return (fig, layers, obsm_layer)

class SCVIPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.scvi", "SCVI", "scvi", order)
        self.dataset = dataset
        self.actions.update(
            fit=FitAction(self.dataset, self.id, loc_class="static"),
        )

    def create_layout(self) -> list:
        cats = self.dataset.get_categoric()
        conts = self.dataset.get_numeric()

        self.components["param-collapse-setup"] = components.CollapsibleDiv(
            page_id=self.id,
            div_id=f"{self.id}-param-collapse-setup",
            collapse_btn_id=f"{self.id}-btn-param-collapse-setup",
            children=[
                html.Div([
                    html.Label("Batch Key"),
                    dcc.Dropdown(
                        options=cats, value=None, id=f"{self.id}-batch_key", clearable=True,
                        placeholder="Select Variable (Optional)",
                        multi=False, style={"flex": "1"},
                    ),
                ], className="param-row-stacked"),
                html.Div([
                    html.Label("Continuous Covariates"),
                    dcc.Dropdown(
                        options=conts, value=None, id=f"{self.id}-continuous_covariate_keys", clearable=True,
                        placeholder="Select Variable(s) (Optional)",
                        multi=True, style={"flex": "1"},
                    ),
                ], className="param-row-stacked"),
                html.Div([
                    html.Label("Categorical Covariates"),
                    dcc.Dropdown(
                        options=cats, value=None, id=f"{self.id}-categorical_covariate_keys", clearable=True,
                        placeholder="Select Variable(s) (Optional)",
                        multi=True, style={"flex": "1"},
                    ),
                ], className="param-row-stacked"),
            ]
        )

        self.components["param-collapse-scvi_model"] = components.CollapsibleDiv(
            page_id=self.id,
            div_id=f"{self.id}-param-collapse-scvi_model",
            collapse_btn_id=f"{self.id}-btn-param-collapse-scvi_model",
            children=components.params_layout(
                SCVI.scvi_model_params, f"{self.id}-model"),
        )
        
        self.components["param-collapse-scvi_train"] = components.CollapsibleDiv(
            page_id=self.id,
            div_id=f"{self.id}-param-collapse-scvi_train",
            collapse_btn_id=f"{self.id}-btn-param-collapse-scvi_train",
            children=components.params_layout(
                SCVI.scvi_train_params, f"{self.id}-train"),
        )

        self.components["param-collapse-scvi_umap"] = components.CollapsibleDiv(
            page_id=self.id,
            div_id=f"{self.id}-param-collapse-scvi_umap",
            collapse_btn_id=f"{self.id}-btn-param-collapse-scvi_umap",
            children=SCVI_UMAP.get_layout(self.id),
        )

        self.components["left_sidebar"] = components.Sidebar(
            page_id=self.id, row="top", side="left",
            title="SCVI Settings", params_children=self._params_layout(),
            apply_btn_id=f"{self.id}-fit-submit", btn_text="Fit SCVI"
        )

        self.components["bot_sidebar"] = components.Sidebar(
            page_id=self.id, row="bot", side="left",
            title="Empty", params_children=[],
            apply_btn_id=None, btn_text="Fit SCVI"
        )

        projection_types = self.dataset.get_scvi_projections()

        type_params = components.FigureHeaderTab(self.id, tab_label="Type", children=[
            html.Div(
                children=[
                    html.Label("Projection Type"),
                    dcc.Dropdown(
                        projection_types, value=next(iter(projection_types), None),
                        id=f"{self.id}-projection-plot-type", clearable=False,
                    ),
                ],
                className="param-row-stacked",
            ),
            # Projection Hue celect
            html.Div(
                children=[
                    html.Label("Color"),
                    dcc.Dropdown(
                        self.dataset.adata.obs_keys()
                        + list(self.dataset.adata.var_names),
                        value=self.dataset.adata.obs_keys()[0],
                        id=f"{self.id}-projection-color",
                        clearable=False,
                    ),
                ],
                className="param-row-stacked",
            ),
        ])

        figure_params = components.FigureHeader(self.id, tabs=[type_params])

        main_figure = html.Div(
            children=[
                html.Div(
                    children= figure_params.create_layout(), className="fig-header",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(id=f"{self.id}-projection-plot",
                                              className="main-plot")
                                )
                            ],
                        )
                    ],
                    id=f"{self.id}-main-body",
                    className="main-body",
                ),
            ],
            className="main",
        )

        secondary_figure = html.Div(
            children=[
                html.Div(
                    children=[
                    ],
                    id=f"{self.id}-secondary-select",
                    className="fig-header",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id=f"{self.id}-secondary-plot",
                                        className="secondary-plot",
                                    )
                                )
                            ],
                        )
                    ],
                    id=f"{self.id}-secondary-body",
                    className="secondary-body",
                ),
            ],
            className="secondary",
        )

        bottom_figure = html.Div(
            children=[],
            id=f"{self.id}-bottom-body",
            className="bottom-body",
        )

        layout = [
            html.Div(
                className="top",
                children=[self.components["left_sidebar"].create_layout(), main_figure, secondary_figure],
            ),
            html.Div(
                className="bottom", children=[self.components["bot_sidebar"].create_layout(), bottom_figure]
            ),
        ]
        return layout

    def _params_layout(self):
        return [
            html.Div([
                dbc.Button(
                    "Setup Parameters", id=f"{self.id}-btn-param-collapse-setup", color="info", n_clicks=0
                ),
                self.components["param-collapse-setup"].create_layout(),
            ], className="param-class"),

            html.Div([
                dbc.Button(
                    "Model Parameters", id=f"{self.id}-btn-param-collapse-scvi_model", color="info", n_clicks=0
                ),
                self.components["param-collapse-scvi_model"].create_layout(),
            ], className="param-class"),

            html.Div([
                dbc.Button(
                    "Train Parameters", id=f"{self.id}-btn-param-collapse-scvi_train", color="info", n_clicks=0
                ),
                self.components["param-collapse-scvi_train"].create_layout(),
            ], className="param-class"),

            html.Div([
                dbc.Button(
                    "Projection Parameters", id=f"{self.id}-btn-param-collapse-scvi_umap", color="info", n_clicks=0
                ),
                self.components["param-collapse-scvi_umap"].create_layout(),
            ], className="param-class"),
        ]
