import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, ctx
from dash.exceptions import PreventUpdate

from ...components import components
from ...components.DashPlot import DashPlot
from ...util.DashAction import DashAction
from . import pca_tools
from ...components.GeneCard import GeneCard
from ...components.CID import CID

import scout

class PlotCorrCircle(DashAction):
    def __init__(
        self, parent_cid: CID, dataset,
        select_pcx_cid: CID,
        select_pcy_cid: CID,
    ):
        super().__init__(parent_cid, dataset)
        self.select_pcx_cid = select_pcx_cid
        self.select_pcy_cid = select_pcy_cid

    def plot(self, pc_x, pc_y):
        fig = scout.ply.pca_corr_circle(
            self.dataset.adata, components=[pc_x, pc_y], layout=pca_tools.default_layout
        )
        return fig

    def setup_callbacks(self, app):
        output = Output(self.parent_cid.to_dict(), "figure")

        # Inputs to Projection
        inputs = dict(
            pc_x=Input(self.select_pcx_cid.to_dict(), "value"),
            pc_y=Input(self.select_pcy_cid.to_dict(), "value"),
            _=Input("url", "pathname")
        )
        @app.dash_app.callback(output=output, inputs=inputs)
        def _(pc_x, pc_y, _):
            return self.plot(pc_x-1, pc_y-1)


class CorrCircle(DashPlot):
    def __init__(self, dataset, page_id, loc_class, select_pcx_cid, select_pcy_cid):
        super().__init__(dataset, page_id, loc_class)
        self.children.update(
            gene_card=GeneCard(self.cid.page_id, self.cid.loc_class, self.dataset)
        )
        self.actions.update(
            plot_correlation_circle=PlotCorrCircle(
                self.cid, self.dataset,
                select_pcx_cid, select_pcy_cid
            )
        )

    def create_layout(self) -> list:
        select_gene_tab = components.FigureHeaderTab(self.page_id, self.loc_class, tab_label="Gene",
            id=f"{self.page_id}-{self.loc_class}-genecard",
            content=self.children["gene_card"].create_layout()
        )

        fig_header = components.FigureHeader(self.page_id, self.loc_class, tabs=[select_gene_tab])

        figure = html.Div([
            html.Div(children=fig_header.create_layout(), className="fig-header"),
            html.Div([
                dcc.Loading(type="circle", children=[
                    html.Div(
                        dcc.Graph(
                            id=self.cid.to_dict(),
                            className=f"{self.loc_class}-plot",
                        )
                    )
                ])
            ], className=f"{self.loc_class}-body"),
        ], className=f"{self.loc_class}")
        return figure