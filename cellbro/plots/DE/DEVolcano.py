from functools import namedtuple

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, ctx
from dash.exceptions import PreventUpdate

from ...components.DashPlot import DashPlot
from ...util.DashAction import DashAction
from ...components import components
from ...components.CID import CID
from ...components.GeneCard import GeneCard
from . import de_tools
from ...components.DropDown import DropDown

import scout

class PlotVolcano(DashAction):
    def __init__(
        self, parent_cid, dataset,
        select_groupby_cid,
        select_reference_cid,
    ):
        super().__init__(parent_cid, dataset)
        self.select_groupby_cid = select_groupby_cid
        self.select_reference_cid = select_reference_cid

    def plot(self, groupby, reference):
        # key = f"rank_genes_{groupby}"

        if not groupby in self.dataset.adata.uns["de"].keys():
            raise PreventUpdate

        if reference is None:
            reference = list(self.dataset.adata.uns["de"][groupby].keys())[0]

        fig = scout.ply.marker_volcano(
            self.dataset.adata.uns["de"][groupby][reference], layout=de_tools.default_layout
        )

        return fig

    def setup_callbacks(self, app):
        outputs = Output(self.parent_cid.to_dict(), "figure")
        inputs = dict(
            groupby=Input(self.select_groupby_cid.to_dict(), "value"),
            reference=Input(self.select_reference_cid.to_dict(), "value"),
        )

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(groupby, reference):

            fig = self.plot(groupby, reference)
            return fig


class DEVolcano(DashPlot):
    def __init__(self, dataset, page_id, loc_class):
        super().__init__(dataset, page_id, loc_class)

        groupby_options = self.dataset.get_rank_genes_groups()
        def get_reference_options():
            groupby_options = self.dataset.get_rank_genes_groups()
            if len(groupby_options) > 0:
                ref_options = de_tools.get_reference_options(self.dataset, groupby_options[0])
            else:
                ref_options = []

            return ref_options

        ref_options = get_reference_options()

        self.children.update(
            gene_card=GeneCard(self.cid.page_id, self.cid.loc_class, self.dataset),
            select_groupby=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-groupby"), 
                options=groupby_options, default=next(iter(groupby_options), None),
                options_callback=lambda: self.dataset.get_rank_genes_groups(),
                update_store_id="de-store"
            ),
            select_reference=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-reference"),
                options=ref_options, default=next(iter(ref_options), None),
                options_callback=lambda: get_reference_options(),
                update_store_id="de-store"
            )
        )
        self.actions.update(
            plot_volcano=PlotVolcano(
                self.cid, self.dataset,
                select_groupby_cid=self.children["select_groupby"].cid,
                select_reference_cid=self.children["select_reference"].cid,
            ),
        )

    def create_layout(self):
        select_ref_tab = components.FigureHeaderTab(
            self.page_id, self.loc_class, tab_label="Reference", content=[
                html.Div([
                    html.Label("Group By"),
                    self.children["select_groupby"].create_layout(),
                    self.children["select_groupby"].get_stores()
                ], className="param-row-stacked"),
                # Projection Hue celect
                html.Div([
                    html.Label("Reference"),
                    self.children["select_reference"].create_layout(),
                    self.children["select_reference"].get_stores()
                ], className="param-row-stacked")
            ]
        )

        select_gene_tab = components.FigureHeaderTab(
            self.page_id, self.loc_class, tab_label="Gene",
            id=f"{self.page_id}-{self.loc_class}-genecard",
            content=self.children["gene_card"].create_layout()
        )

        fig_header = components.FigureHeader(self.page_id, self.loc_class, tabs=[select_ref_tab, select_gene_tab])

        figure = html.Div([
            html.Div(fig_header.create_layout(), className="fig-header"),

            html.Div([
                dcc.Loading(type="circle", children=[
                    html.Div(
                        dcc.Graph(
                            id=self.cid.to_dict(), className=f"{self.loc_class}-plot"
                        )
                    )
                ])
            ], className=f"{self.loc_class}-body"),
        ], className=f"{self.loc_class}")

        return figure
