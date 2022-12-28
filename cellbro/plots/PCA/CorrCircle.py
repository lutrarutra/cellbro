import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

from ...components import components
from ...components.DashFigure import DashFigure
from ...util.DashAction import DashAction
from . import pca_tools
from ...components.GeneListComponents import SelectGene, create_gene_card

import scout


class PlotCorrCircle(DashAction):
    def plot(self, pc_x, pc_y):
        fig = scout.ply.pca_corr_circle(
            self.dataset.adata, components=[pc_x, pc_y], layout=pca_tools.default_layout
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


class CorrCircle(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.actions.update(
            plot_correlation_circle=PlotCorrCircle(self.dataset, self.page_id_prefix),
            select_gene=SelectGene(self.dataset, self.page_id_prefix, self.loc_class),
        )

    def create_layout(self) -> list:
        select_gene_tab = components.FigureHeaderTab(self.page_id_prefix, tab_label="Gene",
            id=f"{self.page_id_prefix}-{self.loc_class}-genecard", children=[
            create_gene_card(None, self.dataset)
        ])

        fig_header = components.FigureHeader(self.page_id_prefix, tabs=[select_gene_tab])

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