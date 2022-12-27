import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

import cellbro.plots.PCA as PCA
import cellbro.util.Components as Components
from cellbro.util.DashPage import DashPage
from cellbro.util.DashAction import DashAction

import scout

from cellbro.plots.DashFigure import DashFigure

class PlotPCA(DashAction):
    def plot(self, color, pc_x, pc_y, continuous_cmap, discrete_cmap):
        fig = scout.ply.projection(
            self.dataset.adata, obsm_layer="X_pca", hue=color, components=[pc_x, pc_y],
            layout=PCA.figure_layout, continuous_cmap=continuous_cmap, discrete_cmap=discrete_cmap
        )
        return fig

    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-main-plot", component_property="figure"),
        ]

        inputs = {
            "color": Input(component_id=f"{self.page_id_prefix}-projection-color", component_property="value"),
            "pc_x": Input(component_id=f"{self.page_id_prefix}-projection-x-component", component_property="value"),
            "pc_y": Input(component_id=f"{self.page_id_prefix}-projection-y-component", component_property="value"),
            "continuous_cmap": Input(
                component_id=f"{self.page_id_prefix}-projection-continuous_cmap", component_property="value"
            ),
            "discrete_cmap": Input(
                component_id=f"{self.page_id_prefix}-projection-discrete_cmap", component_property="value"
            )
        }

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(color, pc_x, pc_y, continuous_cmap, discrete_cmap):
            return [self.plot(color, pc_x-1, pc_y-1, continuous_cmap, discrete_cmap)]

class PlotCorrelationCircle(DashAction):
    def plot(self, pc_x, pc_y):
        fig = scout.ply.pca_correlation_circle(
            self.dataset.adata, components=[pc_x, pc_y], layout=PCA.figure_layout
        )

        return fig

    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-secondary-plot", component_property="figure"),
        ]

        # Inputs to Projection
        inputs = {
            "pc_x": Input(component_id=f"{self.page_id_prefix}-projection-x-component", component_property="value"),
            "pc_y": Input(component_id=f"{self.page_id_prefix}-projection-y-component", component_property="value"),
        }

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(pc_x, pc_y):
            return [self.plot(pc_x-1, pc_y-1)]

class PlotVarianceExplained(DashAction):
    def plot(self, plot_type, n_pcs):
        fig = scout.ply.pca_explain_variance(
            self.dataset.adata, layout=PCA.figure_layout,
            plot_type=plot_type, n_pcs=n_pcs
        )

        return fig

    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-var-plot", component_property="figure")
        ]

        inputs = {
            "plot_type": Input(
                component_id=f"{self.page_id_prefix}-hist-plot_type", component_property="value"
            ),
            "n_pcs": Input(
                component_id=f"{self.page_id_prefix}-hist-n_pcs", component_property="value"
            )
        }

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(plot_type, n_pcs):
            return [self.plot(plot_type, n_pcs)]

class PlotCorrelationExplained(DashAction):
    def plot(self, n_pcs):
        fig = scout.ply.pca_explain_corr(
            self.dataset.adata, layout=PCA.figure_layout,
            n_pcs=n_pcs
        )

        return fig

    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-corr-plot", component_property="figure"),
        ]

        inputs = {
            "n_pcs": Input(
                component_id=f"{self.page_id_prefix}-hist-n_pcs", component_property="value"
            )
        }

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(n_pcs):
            return [self.plot(n_pcs)]

class ExplainVarExplainCorrFigures(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.actions.update(
            plot_var_explained=PlotVarianceExplained(self.dataset, self.page_id_prefix),
            plot_corr_explained=PlotCorrelationExplained(self.dataset, self.page_id_prefix)
        )

    def get_sidebar_params(self) -> list:
        return [
            html.Div(
                children=[
                    html.Label("Plot Type"),
                    dcc.Dropdown(
                        ["Bar", "Line", "Area", "Cumulative"],
                        value="Cumulative",
                        id=f"{self.page_id_prefix}-hist-plot_type",
                        clearable=False,
                    ),
                ],
                className="param-row-stacked",
            ),
            html.Div(
                children=[
                    html.Label("Num. Components"),
                    dcc.Input(
                        id=f"{self.page_id_prefix}-hist-n_pcs", type="number", value=30, min=2, step=1,
                        max=self.dataset.adata.uns["pca"]["variance_ratio"].shape[0] + 1,
                        className="param-input",
                    ),
                ],
                className="param-row",
            ),
        ]

    def create_layout(self) -> list:
        figure = html.Div(
            children=[
                dcc.Loading(
                    type="circle",
                    children=[
                        html.Div(
                            dcc.Graph(id=f"{self.page_id_prefix}-var-plot",className=f"{self.loc_class}-left-plot")
                        )
                    ],
                ),
                dcc.Loading(
                    type="circle",
                    children=[
                        html.Div(
                            dcc.Graph(id=f"{self.page_id_prefix}-corr-plot", className=f"{self.loc_class}-right-plot")
                        )
                    ],
                ),
            ],
            className=f"{self.loc_class}-body",
        )

        return figure

class CorrelationCircleFigure(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.actions.update(
            plot_correlation_circle=PlotCorrelationCircle(self.dataset, self.page_id_prefix),
            select_gene=Components.SelectGene(self.dataset, self.page_id_prefix, self.loc_class),
        )

    def create_layout(self) -> list:
        select_gene_tab = Components.FigureParamTab(self.page_id_prefix, tab_label="Gene", id=f"{self.page_id_prefix}-{self.loc_class}-genecard", children=[
            Components.create_gene_card(self.page_id_prefix, self.loc_class, None, self.dataset)
        ])

        colormap_tab = Components.FigureParamTab(self.page_id_prefix, tab_label="Appearance", children=[
            html.Div([
                html.Label("Continuous Color Map"),
                Components.create_colormap_selector(
                    id=f"{self.page_id_prefix}-projection-continuous_cmap",
                    options=Components.continuous_colormaps,
                    default="viridis",
                )
            ], className="param-row-stacked"),
            html.Div([
                html.Label("Discrete Color Map"),
                Components.create_colormap_selector(
                    id=f"{self.page_id_prefix}-projection-discrete_cmap",
                    options=Components.discrete_colormaps,
                    default="scanpy_default",
                )
            ], className="param-row-stacked")
        ])

        fig_params = Components.FigureParams(self.page_id_prefix, tabs=[select_gene_tab, colormap_tab])

        figure = html.Div(
            children=[
                html.Div(children=fig_params.create_layout(), className="fig-header"),
                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id=f"{self.page_id_prefix}-{self.loc_class}-plot",
                                        className=f"{self.loc_class}-plot",
                                    )
                                )
                            ],
                        )
                    ],
                    className=f"{self.loc_class}-body",
                ),
            ],
            className=f"{self.loc_class}",
        )
        return figure

class PCAProjectionFigure(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.actions.update(
            plot_projection=PlotPCA(self.dataset, self.page_id_prefix),
        )

    def create_layout(self) -> list:
        type_params = Components.FigureParamTab(self.page_id_prefix, tab_label="Type", children=[
            html.Div([
                html.Label("Color"),
                dcc.Dropdown(
                    self.dataset.adata.obs_keys() + list(self.dataset.adata.var_names),
                    value=self.dataset.adata.obs_keys()[0],
                    id=f"{self.page_id_prefix}-projection-color", clearable=False,
                ),
            ], className="param-row-stacked"),
            # X-axis component
            html.Div([
                html.Label("X Component"),
                dbc.Input(
                    id=f"{self.page_id_prefix}-projection-x-component", type="number", value=1, min=1, step=1,
                    max=self.dataset.adata.uns["pca"]["variance_ratio"].shape[0] + 1,
                    className="param-input",
                ),
            ], className="param-row-stacked"),
            # X-axis component
            html.Div(
                children=[
                    html.Label("Y Component"),
                    dbc.Input(
                        id=f"{self.page_id_prefix}-projection-y-component", type="number", value=2, min=1, step=1,
                        max=self.dataset.adata.uns["pca"]["variance_ratio"].shape[0] + 1,
                        className="param-input",
                    ),
                ],
                className="param-row-stacked",
            ),
        ])

        colormap_tab = Components.FigureParamTab(self.page_id_prefix, tab_label="Colormap", children=[
            html.Div([
                html.Label("Continuous Color Map"),
                Components.create_colormap_selector(
                    id=f"{self.page_id_prefix}-projection-continuous_cmap",
                    options=Components.continuous_colormaps,
                    default="viridis",
                )
            ], className="param-row-stacked"),
            html.Div([
                html.Label("Discrete Color Map"),
                Components.create_colormap_selector(
                    id=f"{self.page_id_prefix}-projection-discrete_cmap",
                    options=Components.discrete_colormaps,
                    default="scanpy_default",
                )
            ], className="param-row-stacked"),

        ])

        figure_params = Components.FigureParams(self.page_id_prefix, tabs=[type_params, colormap_tab])

        figure = html.Div(
            children=[
                html.Div(
                    children=figure_params.create_layout(),
                    className="fig-header",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id=f"{self.page_id_prefix}-{self.loc_class}-plot",
                                        className=f"{self.loc_class}-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    className=f"{self.loc_class}-body",
                ),
            ],
            className="main",
        )
        return figure

class PCAPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.pca", "PCA", "pca", order)
        self.dataset = dataset
        self.components.update(
            pca_figure=PCAProjectionFigure(self.dataset, self.id, "main"),
            correlation_circle_figure=CorrelationCircleFigure(self.dataset, self.id, "secondary"),
            bottom_figure=ExplainVarExplainCorrFigures(self.dataset, self.id, "bottom")
        )

    def create_layout(self) -> list:
        self.components["left_sidebar"] = Components.Sidebar(
            page_id_prefix=self.id, row="top", side="left",
            title="PCA Projection Settings",
            params_children=self.components["pca_figure"].get_sidebar_params(),
            apply_btn_id=None, btn_text="Plot"
        )

        self.components["bot_sidebar"] = Components.Sidebar(
            page_id_prefix=self.id, row="bot", side="left",
            title="PCA Plots",
            params_children=self.components["bottom_figure"].get_sidebar_params(),
            apply_btn_id=None, btn_text="Plot"
        )

        main_figure = self.components["pca_figure"].create_layout()
        secondary_figure = self.components["correlation_circle_figure"].create_layout()
        bottom_figure = self.components["bottom_figure"].create_layout()

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

    # def _top_params_layout(self):
        # divs = []
        # for key, param in PCA.pca_params.items():
        #     divs.append(
        #         html.Div(
        #             children=[
        #                 html.Label(param.name, className="param-label",),
        #                 dcc.Input(
        #                     id=f"{self.id}-{key}", type=param.input_type, value=param.value,
        #                     step=param.step if param.step != None else 0.1,
        #                     className="param-input",
        #                 ),
        #             ],
        #             className="param-row",
        #         )
        #     )

        # # layout = html.Div(children=[divs])
        # return divs
