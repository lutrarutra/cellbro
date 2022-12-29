from dash import Input, Output, State, ctx, dcc, html
from dash.exceptions import PreventUpdate

from ...components.DashFigure import DashFigure
from ...components import components
from ...util.DashAction import DashAction
from ...components.GeneListComponents import SelectGene, create_gene_card, UpdateGeneList
from .qc_tools import default_layout

import scout

class Plot(DashAction):
    def plot(self):
        fig = scout.ply.dispersion_plot(
            self.dataset.adata, layout=default_layout
        )
        return fig

    def setup_callbacks(self, app):
        output = Output(f"{self.page_id_prefix}-{self.loc_class}-plot", "figure")
        inputs = dict(
            submit=Input(f"{self.page_id_prefix}-apply-btn", "n_clicks"),
        )
        @app.dash_app.callback(output=output, inputs=inputs)
        def _(submit):
            if self.dataset.qc_done():
                return self.plot()

            raise PreventUpdate


class DispersionPlot(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.actions.update(
            select_gene=SelectGene(dataset, self.page_id_prefix, self.loc_class),
            update_gene_list=UpdateGeneList(self.dataset, self.page_id_prefix, self.loc_class),
            plot=Plot(dataset, page_id_prefix, loc_class)
        )

    def create_layout(self):

        select_gene_tab = components.FigureHeaderTab(
            self.page_id_prefix, tab_label="Gene", id=f"{self.page_id_prefix}-{self.loc_class}-genecard",
            children=create_gene_card(self.page_id_prefix, self.loc_class, None, self.dataset)
        )

        figure_header = components.FigureHeader(
            self.page_id_prefix, [select_gene_tab]
        )

        figure = html.Div([
            html.Div(children=figure_header.create_layout(), className="fig-header"),
            html.Div([
                dcc.Loading(type="circle", children=[
                    html.Div(
                        dcc.Graph(id=f"{self.page_id_prefix}-{self.loc_class}-plot", className=f"{self.loc_class}-plot")
                    )
                ])
            ], className=f"{self.loc_class}-body")
        ], className=f"{self.loc_class}")

        return figure
