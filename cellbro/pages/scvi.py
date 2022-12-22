import dash
from dash import Input, Output, State, dcc, html, ctx
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

import cellbro.plots.SCVI as scvi_plots
import cellbro.util.Components as Components
import cellbro.plots.Projection as Projection
from cellbro.util.DashPage import DashPage
from cellbro.util.DashAction import DashAction
from cellbro.plots.UMAP import SCVI_UMAP

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
            Output(f"{self.page_id_prefix}-projection-plot", "figure"),
            Output(f"{self.page_id_prefix}-projection-plot-type", "options"),
            Output(component_id=f"{self.page_id_prefix}-projection-plot-type", component_property="value"),
        ]
        inputs = {
            "submit": Input(f"{self.page_id_prefix}-fit-submit", "n_clicks"),
            "color": Input(
                component_id=f"{self.page_id_prefix}-projection-color", component_property="value"
            ),
            "obsm_layer": Input(
                component_id=f"{self.page_id_prefix}-projection-plot-type", component_property="value"
            ),
        }
        state = {
            "setup_continuous_covariate_keys": Input(f"{self.page_id_prefix}-continuous_covariate_keys", "value"),
            "setup_categorical_covariate_keys": Input(f"{self.page_id_prefix}-categorical_covariate_keys", "value"),
            "setup_batch_key": Input(f"{self.page_id_prefix}-batch_key", "value"),
        }
        for key, param in scvi_plots.scvi_model_params.items():
            state[f"model_{param.key}"] = State(f"{self.page_id_prefix}-model-{param.key}", "value")

        for key, param in scvi_plots.scvi_train_params.items():
            state[f"train_{param.key}"] = State(f"{self.page_id_prefix}-train-{param.key}", "value")

        for key in SCVI_UMAP._params.keys():
            state[f"scvi_umap_{key}"] = State(
                component_id=f"{self.page_id_prefix}-projection-scvi_umap-{key}", component_property="value"
            )

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, color, obsm_layer, **kwargs):
            if ctx.triggered_id == f"{self.page_id_prefix}-fit-submit":
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
            fit=FitAction(self.dataset, self.id),
        )

    def create_layout(self) -> list:
        cats = self.dataset.get_categoric()
        conts = self.dataset.get_numeric()

        self.actions["param-collapse-setup"] = Components.CollapseDiv(
            page_id_prefix=self.id,
            id=f"{self.id}-param-collapse-setup",
            btn_id=f"{self.id}-btn-param-collapse-setup",
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

        self.actions["param-collapse-scvi_model"] = Components.CollapseDiv(
            page_id_prefix=self.id,
            id=f"{self.id}-param-collapse-scvi_model",
            btn_id=f"{self.id}-btn-param-collapse-scvi_model",
            children=Components.params_layout(
                scvi_plots.scvi_model_params, f"{self.id}-model"),
        )
        
        self.actions["param-collapse-scvi_train"] = Components.CollapseDiv(
            page_id_prefix=self.id,
            id=f"{self.id}-param-collapse-scvi_train",
            btn_id=f"{self.id}-btn-param-collapse-scvi_train",
            children=Components.params_layout(
                scvi_plots.scvi_train_params, f"{self.id}-train"),
        )

        self.actions["param-collapse-scvi_umap"] = Components.CollapseDiv(
            page_id_prefix=self.id,
            id=f"{self.id}-param-collapse-scvi_umap",
            btn_id=f"{self.id}-btn-param-collapse-scvi_umap",
            children=Projection.Projection.get_layout(SCVI_UMAP, self.id),
        )

        top_sidebar = Components.create_sidebar(
            id=f"{self.id}-top-sidebar", class_name="top-sidebar",
            title="SCVI Settings",
            params_children=self._params_layout(),
            btn_id=f"{self.id}-fit-submit", btn_text="Fit SCVI"
        )

        bot_sidebar = Components.create_sidebar(
            id=f"{self.id}-bot-sidebar", class_name="bot-sidebar",
            title="Empty",
            params_children=[],
            btn_id=None, btn_text="Fit SCVI"
        )

        projection_types = self.dataset.get_scvi_projections()

        main_figure = html.Div(
            children=[
                html.Div(
                    children=[
                        html.Div(
                            children=[
                                html.Label("Projection Type"),
                                dcc.Dropdown(
                                    projection_types, value=next(iter(projection_types), None),
                                    id=f"{self.id}-projection-plot-type", clearable=False,
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
                                    id=f"{self.id}-projection-color",
                                    clearable=False,
                                ),
                            ],
                            className="param-column",
                        ),
                    ],
                    id=f"{self.id}-main-select",
                    className="top-parameters",
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
                    id=f"{self.id}-main-figure",
                    className="main-figure",
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
                    className="top-parameters",
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
                    id=f"{self.id}-secondary-figure",
                    className="secondary-figure",
                ),
            ],
            className="secondary",
        )

        bottom_figure = html.Div(
            children=[],
            id=f"{self.id}-bottom-figure",
            className="bottom-figure",
        )

        layout = [
            html.Div(
                className="top",
                children=[top_sidebar, main_figure, secondary_figure],
            ),
            html.Div(
                className="bottom", children=[bot_sidebar, bottom_figure]
            ),
        ]
        return layout

    def _params_layout(self):
        return [
            html.Div([
                dbc.Button(
                    "Setup Parameters", id=f"{self.id}-btn-param-collapse-setup", color="info", n_clicks=0
                ),
                self.actions["param-collapse-setup"].layout,
            ], className="param-class"),

            html.Div([
                dbc.Button(
                    "Model Parameters", id=f"{self.id}-btn-param-collapse-scvi_model", color="info", n_clicks=0
                ),
                self.actions["param-collapse-scvi_model"].layout,
            ], className="param-class"),

            html.Div([
                dbc.Button(
                    "Train Parameters", id=f"{self.id}-btn-param-collapse-scvi_train", color="info", n_clicks=0
                ),
                self.actions["param-collapse-scvi_train"].layout,
            ], className="param-class"),

            html.Div([
                dbc.Button(
                    "Projection Parameters", id=f"{self.id}-btn-param-collapse-scvi_umap", color="info", n_clicks=0
                ),
                self.actions["param-collapse-scvi_umap"].layout,
            ], className="param-class"),
        ]
