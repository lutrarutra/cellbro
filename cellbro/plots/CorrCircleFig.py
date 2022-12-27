import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

from ..util import Components
from ..util.DashFigure import DashFigure
from ..util.DashAction import DashAction
from . import PCA

import scout


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


class CorrCircleFig(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.actions.update(
            plot_correlation_circle=PlotCorrelationCircle(self.dataset, self.page_id_prefix),
            select_gene=Components.SelectGene(self.dataset, self.page_id_prefix, self.loc_class),
        )

    def create_layout(self) -> list:
        select_gene_tab = Components.FigureParamTab(self.page_id_prefix, tab_label="Gene",
            id=f"{self.page_id_prefix}-{self.loc_class}-genecard", children=[
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

        fig_header = Components.FigureParams(self.page_id_prefix, tabs=[select_gene_tab, colormap_tab])

        figure = html.Div(
            children=[
                html.Div(children=fig_header.create_layout(), className="fig-header"),
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