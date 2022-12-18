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

        return self.plot(params)

    def plot(self, params: dict):
        if "neighbors_scvi" not in self.dataset.get_neighbors():
            raise PreventUpdate
        
        _, kwargs = Projection.parse_params(params)
        projection = SCVI_UMAP(dataset=self.dataset, **kwargs)
        projection.apply()

        return [projection.plot()]

    def setup_callbacks(self, app):
        output = [
            Output("scvi-projection-plot", "figure"),
        ]
        inputs = {
            "submit": Input("fit-submit", "n_clicks"),
            "projection_color": Input(
                component_id="scvi-projection-color", component_property="value"
            ),
            "projection_type": Input(
                component_id="scvi-projection-type", component_property="value"
            ),
        }
        state = {
            "setup_continuous_covariate_keys": Input("continuous_covariate_keys", "value"),
            "setup_categorical_covariate_keys": Input("categorical_covariate_keys", "value"),
            "setup_batch_key": Input("batch_key", "value"),
        }
        for key, param in scvi_plots.scvi_model_params.items():
            state[f"model_{param.key}"] = State(f"scvi-model-{param.key}", "value")

        for key, param in scvi_plots.scvi_train_params.items():
            state[f"train_{param.key}"] = State(f"scvi-train-{param.key}", "value")

        for key in SCVI_UMAP._params.keys():
            state[f"scvi_umap_{key}"] = State(
                component_id=f"projection-scvi_umap-{key}", component_property="value"
            )

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(**kwargs):
            _ = kwargs.pop("submit")
            if ctx.triggered_id == "fit-submit":
                return self.apply(params=kwargs)
            # if submit is None:
            #     raise PreventUpdate
            
            return self.plot(params=kwargs)

class SCVIPage(DashPage):
    def __init__(self, dataset, app, order):
        super().__init__("pages.scvi", "SCVI", "/scvi", order)
        self.dataset = dataset
        self.actions.update(
            fit=FitAction(self.dataset),
        )
        self.layout = self.create_layout()
        self.setup_callbacks(app)


    def create_layout(self) -> list:
        cats = self.dataset.get_categoric()
        conts = self.dataset.get_numeric()

        self.actions["param-collapse-scvi-setup"] = Components.CollapseDiv(
            id="param-collapse-scvi-setup",
            btn_id="btn-param-collapse-scvi-setup",
            children=[
                html.Div([
                    html.Label("Batch Key"),
                    dcc.Dropdown(
                        options=cats, value=None, id="batch_key", clearable=True,
                        placeholder="Select Variable (Optional)",
                        multi=False, style={"flex": "1"},
                    ),
                ], className="param-row-stacked"),
                html.Div([
                    html.Label("Continuous Covariates"),
                    dcc.Dropdown(
                        options=conts, value=None, id="continuous_covariate_keys", clearable=True,
                        placeholder="Select Variable(s) (Optional)",
                        multi=True, style={"flex": "1"},
                    ),
                ], className="param-row-stacked"),
                html.Div([
                    html.Label("Categorical Covariates"),
                    dcc.Dropdown(
                        options=cats, value=None, id="categorical_covariate_keys", clearable=True,
                        placeholder="Select Variable(s) (Optional)",
                        multi=True, style={"flex": "1"},
                    ),
                ], className="param-row-stacked"),
            ]
        )

        self.actions["param-collapse-scvi_model"] = Components.CollapseDiv(
            id="param-collapse-scvi_model",
            btn_id="btn-param-collapse-scvi_model",
            children=Components.params_layout(
                scvi_plots.scvi_model_params, "scvi-model"),
        )
        
        self.actions["param-collapse-scvi_train"] = Components.CollapseDiv(
            id="param-collapse-scvi_train",
            btn_id="btn-param-collapse-scvi_train",
            children=Components.params_layout(
                scvi_plots.scvi_train_params, "scvi-train"),
        )

        self.actions["param-collapse-scvi_umap"] = Components.CollapseDiv(
            id="param-collapse-scvi_umap",
            btn_id="btn-param-collapse-scvi_umap",
            children=Projection.Projection.get_layout(SCVI_UMAP),
        )

        top_sidebar = Components.create_sidebar(
            id="scvi-top-sidebar", class_name="top-sidebar",
            title="SCVI Settings",
            params_children=self._params_layout(),
            btn_id="fit-submit", btn_text="Rune"
        )

        bot_sidebar = Components.create_sidebar(
            id="scvi-bot-sidebar", class_name="bot-sidebar",
            title="Empty",
            params_children=[],
            btn_id=None, btn_text="Rune"
        )

        main_figure = html.Div(
            children=[
                html.Div(
                    children=[
                        html.Div(
                            children=[
                                html.Label("Projection Type"),
                                dcc.Dropdown(
                                    ["SCVI-UMAP"], value="SCVI-UMAP",
                                    id="scvi-projection-type", clearable=False,
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
                                    id="scvi-projection-color",
                                    clearable=False,
                                ),
                            ],
                            className="param-column",
                        ),
                    ],
                    id="scvi-main-select",
                    className="top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            id="loading-scvi-main",
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(id="scvi-projection-plot",
                                              className="main-plot")
                                )
                            ],
                        )
                    ],
                    id="scvi-main-figure",
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
                    id="scvi-secondary-select",
                    className="top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            id="loading-scvi-secondary",
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id="scvi-secondary-plot",
                                        className="secondary-plot",
                                    )
                                )
                            ],
                        )
                    ],
                    id="scvi-secondary-figure",
                    className="secondary-figure",
                ),
            ],
            className="secondary",
        )

        bottom_figure = html.Div(
            children=[],
            id="scvi-bottom-figure",
            className="bottom-figure",
        )

        layout = [
            html.Div(
                id="top",
                className="top",
                children=[top_sidebar, main_figure, secondary_figure],
            ),
            html.Div(
                id="bottom", className="bottom", children=[bot_sidebar, bottom_figure]
            ),
        ]
        return layout

    def _params_layout(self):
        return [
            html.Div([
                dbc.Button(
                    "Setup Parameters", id="btn-param-collapse-scvi-setup", color="info", n_clicks=0
                ),
                self.actions["param-collapse-scvi-setup"].layout,
            ], className="param-class"),

            html.Div([
                dbc.Button(
                    "Model Parameters", id="btn-param-collapse-scvi_model", color="info", n_clicks=0
                ),
                self.actions["param-collapse-scvi_model"].layout,
            ], className="param-class"),

            html.Div([
                dbc.Button(
                    "Train Parameters", id="btn-param-collapse-scvi_train", color="info", n_clicks=0
                ),
                self.actions["param-collapse-scvi_train"].layout,
            ], className="param-class"),

            html.Div([
                dbc.Button(
                    "Projection Parameters", id="btn-param-collapse-scvi_umap", color="info", n_clicks=0
                ),
                self.actions["param-collapse-scvi_umap"].layout,
            ], className="param-class"),
        ]
