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

    def setup_callbacks(self, app):
        outputs = [Output(f"{self.page_id_prefix}-secondary-select", "children")]
        inputs = {
            "click_data": Input(f"{self.page_id_prefix}-secondary-plot", "clickData"),
        }

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(**kwargs):
            if kwargs["click_data"] is None:
                raise PreventUpdate
            return self.apply(params=kwargs)

class PlotAction(DashAction):
    def apply(self, params):
        return PCA.PCA(self.dataset, **params).plot()

    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-main-plot", component_property="figure"),
            Output(component_id=f"{self.page_id_prefix}-secondary-plot", component_property="figure"),
            Output(component_id=f"{self.page_id_prefix}-var-plot", component_property="figure"),
            Output(component_id=f"{self.page_id_prefix}-corr-plot", component_property="figure"),
        ]

        # Inputs to Projection
        inputs = {
            "color": Input(component_id=f"{self.page_id_prefix}-projection-color", component_property="value"),
            "pc_x": Input(component_id=f"{self.page_id_prefix}-projection-x-component", component_property="value"),
            "pc_y": Input(component_id=f"{self.page_id_prefix}-projection-y-component", component_property="value"),
            "hist_type": Input(
                component_id=f"{self.page_id_prefix}-hist_type", component_property="value"
            ),
            "hist_n_pcs": Input(
                component_id=f"{self.page_id_prefix}-hist_n_pcs", component_property="value"
            ),
        }

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(**kwargs):
            return self.apply(params=kwargs)


class PCAPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.pca", "PCA", "pca", order)
        self.dataset = dataset
        self.actions.update(
            plot_action=PlotAction(self.dataset, page_id_prefix=self.id),
            click_action=ClickAction(self.dataset, page_id_prefix=self.id)
        )

    def create_layout(self) -> list:
        top_sidebar = Components.create_sidebar(
            id=f"{self.id}-top-sidebar", class_name="top-sidebar",
            title="PCA Projection Settings",
            params_children=self._top_params_layout(),
            btn_id=None, btn_text="Plot"
        )

        bot_sidebar = Components.create_sidebar(
            id=f"{self.id}-bot-sidebar", class_name="bot-sidebar",
            title="PCA Plots",
            params_children=self._bot_params_layout(),
            btn_id=None, btn_text="Plot"
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
                                    self.dataset.adata.obs_keys() + list(self.dataset.adata.var_names),
                                    value=self.dataset.adata.obs_keys()[0],
                                    id=f"{self.id}-projection-color", clearable=False,
                                ),
                            ],
                            className="param-column",
                        ),
                        # X-axis component
                        html.Div(
                            children=[
                                html.Label("X Component"),
                                dcc.Input(
                                    id=f"{self.id}-projection-x-component", type="number", value=1, min=1, step=1,
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
                                    id=f"{self.id}-projection-y-component", type="number", value=2, min=1, step=1,
                                    max=self.dataset.adata.uns["pca"]["variance_ratio"].shape[0] + 1,
                                    className="param-input",
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
                                    dcc.Graph(id=f"{self.id}-main-plot",
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
                        Components.create_gene_card(None, self.dataset),
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
            children=[
                dcc.Loading(
                    type="circle",
                    children=[
                        html.Div(
                            dcc.Graph(id=f"{self.id}-var-plot",className="bottom-left-plot")
                        )
                    ],
                ),
                dcc.Loading(
                    type="circle",
                    children=[
                        html.Div(
                            dcc.Graph(id=f"{self.id}-corr-plot", className="bottom-right-plot")
                        )
                    ],
                ),
            ],
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

    def _bot_params_layout(self):
        return [
            html.Div(
                children=[
                    html.Label("Plot Type"),
                    dcc.Dropdown(
                        ["Bar", "Linefig", "Area", "Cumulative"],
                        value="Cumulative",
                        id=f"{self.id}-hist_type",
                        clearable=False,
                    ),
                ],
                className="param-row-stacked",
            ),
            html.Div(
                children=[
                    html.Label("Num. Components"),
                    dcc.Input(
                        id=f"{self.id}-hist_n_pcs", type="number", value=30, min=2, step=1,
                        max=self.dataset.adata.uns["pca"]["variance_ratio"].shape[0] + 1,
                        className="param-input",
                    ),
                ],
                className="param-row",
            ),
        ]

    def _top_params_layout(self):
        divs = []
        for key, param in PCA.pca_params.items():
            divs.append(
                html.Div(
                    children=[
                        html.Label(param.name, className="param-label",),
                        dcc.Input(
                            id=f"{self.id}-{key}", type=param.input_type, value=param.value,
                            step=param.step if param.step != None else 0.1,
                            className="param-input",
                        ),
                    ],
                    className="param-row",
                )
            )

        # layout = html.Div(children=[divs])
        return divs
