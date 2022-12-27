import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

import cellbro.plots.PCA as PCA
import cellbro.util.Components as Components
from cellbro.util.DashPage import DashPage
from cellbro.util.DashAction import DashAction

from ..plots import projection as prj
from ..plots.CorrCircleFig import CorrCircleFig

import scout

from cellbro.util.DashFigure import DashFigure

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

class PCAPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.pca", "PCA", "pca", order)
        self.dataset = dataset
        self.components.update(
            pca_figure=prj.PCAProjectionFigure(self.dataset, self.id, "main"),
            correlation_circle_figure=CorrCircleFig(self.dataset, self.id, "secondary"),
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
