from dash import Input, Output, State, ctx, dcc, html
from dash.exceptions import PreventUpdate

from ...components.DashPlot import DashPlot
from ...components import components
from ...util.DashAction import DashAction
from .qc_tools import default_layout
from ...components.GeneCard import GeneCard

import scout

class Plot(DashAction):
    def plot(self):
        fig = scout.ply.dispersion_plot(
            self.dataset.adata, layout=default_layout
        )
        return fig

    def setup_callbacks(self, app):
        output = Output(self.parent_cid.to_str(), "figure")
        inputs = dict(
            submit=Input(f"{self.page_id}-main-sidebar-apply_btn", "n_clicks"),
        )
        @app.dash_app.callback(output=output, inputs=inputs)
        def _(submit):
            if self.dataset.qc_done():
                return self.plot()

            raise PreventUpdate


class DispersionPlot(DashPlot):
    def __init__(self, dataset, page_id, loc_class):
        super().__init__(dataset, page_id, loc_class)
        self.actions.update(
            # select_gene=SelectGene(self.cid, dataset),
            # update_genelist=UpdateGeneList(self.cid, self.dataset),
            plot=Plot(self.cid, dataset)
        )
        self.children.update(
            gene_card=GeneCard(self.cid.page_id, self.cid.loc_class, self.dataset)
        )

    def create_layout(self):
        select_gene_tab = components.FigureHeaderTab(
            self.page_id, self.loc_class, tab_label="Gene", id=f"{self.page_id}-{self.loc_class}-genecard",
            content=self.children["gene_card"].create_layout()
        )

        figure_header = components.FigureHeader(
            self.page_id, self.loc_class, [select_gene_tab]
        )

        figure = html.Div([
            html.Div(children=figure_header.create_layout(), className="fig-header"),
            html.Div([
                dcc.Loading(type="circle", children=[
                    html.Div(
                        dcc.Graph(id=f"{self.page_id}-{self.loc_class}-plot", className=f"{self.loc_class}-plot")
                    )
                ])
            ], className=f"{self.loc_class}-body")
        ], className=f"{self.loc_class}")

        return figure
