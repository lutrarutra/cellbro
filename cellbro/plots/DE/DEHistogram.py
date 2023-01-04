from functools import namedtuple

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, ctx
from dash.exceptions import PreventUpdate

from ...components.DashPlot import DashPlot
from ...util.DashAction import DashAction
from ...components.CID import CID
from . import de_tools
from ...components.DropDown import DropDown
from . import de_tools 
from ...components import components

import scout

class Plot(DashAction):
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
        fig = scout.ply.pval_histogram(
            self.dataset.adata.uns["de"][groupby][reference], layout=de_tools.default_layout
        )
        return fig

    def setup_callbacks(self, app):
        outputs = [
            Output(self.parent_cid.to_dict(), "figure"),
        ]
        inputs = dict(
            groupby=Input(self.select_groupby_cid.to_dict(), "value"),
            reference=Input(self.select_reference_cid.to_dict(), "value"),
        )

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(groupby, reference):
            if groupby is None or reference is None:
                raise PreventUpdate

            return self.plot(groupby, reference)


class DEHistogram(DashPlot):
    def __init__(self, dataset, page_id, loc_class, select_groupby_cid=None, select_reference_cid=None):
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
        
        if select_groupby_cid is None:
            assert False, "Not implemented"
            self.select_reference_cid = CID(self.page_id, self.loc_class, "select-groupby")
            self.children.update(
                select_groupby=DropDown(
                    cid=self.select_reference_cid,
                    options=groupby_options, default=next(iter(groupby_options), None),
                    options_callback=lambda: self.dataset.get_rank_genes_groups(),
                    update_store_id="de-store"
                )
            )
        else:
            self.select_groupby_cid = select_groupby_cid

        if select_reference_cid is None:
            assert False, "Not implemented"
            self.select_groupby_cid = CID(self.page_id, self.loc_class, "select-groupby")
            self.children.update(
                select_reference=DropDown(
                    cid=self.select_groupby_cid,
                    options=ref_options, default=next(iter(ref_options), None),
                    options_callback=lambda: get_reference_options(),
                    update_store_id="de-store"
                )
            )
        else:
            self.select_reference_cid = select_reference_cid


        self.children.update(
            select_plot_type=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-plot_type"),
                options=["pval_histogram"], default="pval_histogram",
            )
        )
        
        self.actions.update(
            plot=Plot(
                parent_cid=self.cid, dataset=self.dataset,
                select_groupby_cid=self.select_groupby_cid,
                select_reference_cid=self.select_reference_cid,
            )
        )

    def create_layout(self):
        plot_type_params = components.FigureHeaderTab(self.page_id, self.loc_class, tab_label="Type", content=[
            # Volcano Group By Select
            html.Div([
                html.Label("Plot Type"),
                self.children["select_plot_type"].create_layout(),
                self.children["select_plot_type"].get_stores(),
            ], className="param-row-stacked")
        ])

        figure_params = components.FigureHeader(self.page_id, self.loc_class, tabs=[plot_type_params])

        figure = html.Div([
            html.Div(figure_params.create_layout(), className="fig-header"),
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
        ], className=self.loc_class.name)
        return figure
        
        

