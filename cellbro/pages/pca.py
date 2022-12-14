import dash
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

import cellbro.plots.PCA as PCA
import cellbro.util.Components as Components
from cellbro.util.DashPage import DashPage
from cellbro.util.DashAction import DashAction

# TODO: reusable class
class ClickAction(DashAction):
    def apply(self, params):
        gene = params["click_data"]["points"][0]["hovertext"]
        element = Components.create_gene_card(gene, self.dataset)
        return [element]

    def setup_callbacks(self, dash_app):
        outputs = [Output("pca-secondary-select", "children")]
        inputs = {
            "click_data": Input("pca-secondary-plot", "clickData"),
        }

        @dash_app.callback(output=outputs, inputs=inputs)
        def _(**kwargs):
            if kwargs["click_data"] is None:
                raise PreventUpdate
            return self.apply(params=kwargs)

class PlotAction(DashAction):
    def apply(self, params):
        return PCA.PCA(self.dataset, **params).plot()

    def setup_callbacks(self, dash_app):
        output = [
            Output(component_id="pca-main-plot", component_property="figure"),
            Output(component_id="pca-secondary-plot", component_property="figure"),
            Output(component_id="pca-var-plot", component_property="figure"),
            Output(component_id="pca-corr-plot", component_property="figure"),
        ]

        # Inputs to Projection
        inputs = {
            "color": Input(component_id="pca-color", component_property="value"),
            "pc_x": Input(component_id="pca-x-component", component_property="value"),
            "pc_y": Input(component_id="pca-y-component", component_property="value"),
            "hist_type": Input(
                component_id="pca-hist_type", component_property="value"
            ),
            "hist_n_pcs": Input(
                component_id="pca-hist_n_pcs", component_property="value"
            ),
        }

        @dash_app.callback(output=output, inputs=inputs)
        def _(**kwargs):
            return self.apply(params=kwargs)


class PCAPage(DashPage):
    def __init__(self, dataset, dash_app, order):
        super().__init__("pages.pca", "PCA", "/pca", order)
        self.dataset = dataset
        self.layout = self.create_layout()
        self.actions = dict(
            plot_action=PlotAction(self.dataset),
            click_action=ClickAction(self.dataset)
        )
        self.setup_callbacks(dash_app)

    def create_layout(self) -> list:
        top_sidebar = html.Div(
            children=[
                html.Div(
                    [
                        html.H3("PCA Projection Settings"),
                    ],
                    className="sidebar-header",
                ),
                dcc.Loading(
                    type="circle",
                    children=[
                        html.Div(children=[self._params_layout()], className="sidebar-parameters"),
                        html.Div([
                                # dbc.Button("Filter", color="primary", className="mr-1", id="pca-submit"),
                        ], className="sidebar-footer",),
                    ],
                ),
            ],
            className="top-sidebar sidebar",
        )

        bottom_sidebar = html.Div(
            children=[
                html.Div([
                    html.H3("PCA Plots"),
                ], className="sidebar-header"),
                dcc.Loading(
                    type="circle",
                    children=[
                        # Num components
                        html.Div(
                            children=[
                                html.Div(
                                    children=[
                                        html.Label("Plot Type"),
                                        dcc.Dropdown(
                                            ["Bar", "Linefig", "Area", "Cumulative"],
                                            value="Bar",
                                            id="pca-hist_type",
                                            clearable=False,
                                        ),
                                    ],
                                    className="param-row-stacked",
                                ),
                                html.Div(
                                    children=[
                                        html.Label("Num. Components"),
                                        dcc.Input(
                                            id="pca-hist_n_pcs", type="number", value=30, min=2, step=1,
                                            max=self.dataset.adata.uns["pca"]["variance_ratio"].shape[0] + 1,
                                            className="param-input",
                                        ),
                                    ],
                                    className="param-row",
                                ),
                            ],
                            className="sidebar-parameters",
                        ),
                    ],
                ),
            ],
            id="pca-bottom-sidebar",
            className="bottom-sidebar sidebar",
        )

        main_figure = html.Div(
            children=[
                html.Div(
                    children=[
                        # Projection Color
                        html.Div(
                            children=[
                                html.Label("Color"),
                                dcc.Dropdown(
                                    self.dataset.adata.obs_keys() + self.dataset.adata.var_names.tolist(),
                                    value=self.dataset.adata.obs_keys()[0],
                                    id="pca-color", clearable=False,
                                ),
                            ],
                            className="param-column",
                        ),
                        # X-axis component
                        html.Div(
                            children=[
                                html.Label("X Component"),
                                dcc.Input(
                                    id="pca-x-component", type="number", value=1, min=1, step=1,
                                    max=self.dataset.adata.uns["pca"]["variance_ratio"].shape[0] + 1,
                                    className="param-input",
                                ),
                            ],
                            className="param-column",
                        ),
                        # X-axis component
                        html.Div(
                            children=[
                                html.Label("Y Component"),
                                dcc.Input(
                                    id="pca-y-component", type="number", value=2, min=1, step=1,
                                    max=self.dataset.adata.uns["pca"]["variance_ratio"].shape[0] + 1,
                                    className="param-input",
                                ),
                            ],
                            className="param-column",
                        ),
                    ],
                    id="pca-main-select",
                    className="top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            id="loading-pca-main",
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(id="pca-main-plot",
                                              className="main-plot")
                                )
                            ],
                        )
                    ],
                    id="pca-main-figure",
                    className="main-figure",
                ),
            ],
            className="main",
        )

        secondary_figure = html.Div(
            children=[
                html.Div(
                    children=[
                        Components.create_gene_card(None, self.dataset),
                    ],
                    id="pca-secondary-select",
                    className="top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            id="loading-pca-secondary",
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id="pca-secondary-plot",
                                        className="secondary-plot",
                                    )
                                )
                            ],
                        )
                    ],
                    id="pca-secondary-figure",
                    className="secondary-figure",
                ),
            ],
            className="secondary",
        )

        bottom_figure = html.Div(
            children=[
                dcc.Loading(
                    id="loading-var-bottom",
                    className="loading-bottom",
                    type="circle",
                    children=[
                        html.Div(
                            dcc.Graph(id="pca-var-plot",
                                      className="bottom-left-plot")
                        )
                    ],
                ),
                dcc.Loading(
                    id="loading-corr-bottom",
                    className="loading-bottom",
                    type="circle",
                    children=[
                        html.Div(
                            dcc.Graph(id="pca-corr-plot",
                                      className="bottom-right-plot")
                        )
                    ],
                ),
            ],
            id="pca-bottom-figure",
            className="bottom-figure",
        )

        layout = [
            html.Div(
                id="top",
                className="top",
                children=[top_sidebar, main_figure, secondary_figure],
            ),
            html.Div(
                id="bottom", className="bottom", children=[bottom_sidebar, bottom_figure]
            ),
        ]
        return layout

    def _params_layout(self):
        divs = []
        for key, param in PCA.pca_params.items():
            divs.append(
                html.Div(
                    children=[
                        html.Label(param.name, className="param-label",),
                        dcc.Input(
                            id=f"pca-{key}", type=param.input_type, value=param.value,
                            step=param.step if param.step != None else 0.1,
                            className="param-input",
                        ),
                    ],
                    className="param-row",
                )
            )

        layout = html.Div(children=[divs])
        return layout
