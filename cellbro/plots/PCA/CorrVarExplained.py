import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

from ...components import components
from ...components.DashPlot import DashPlot
from ...util.DashAction import DashAction
from . import pca_tools
from ...components.GeneCard import GeneCard
from ...components.CID import CID
from ...components.InputField import InputField
from .pca_tools import default_layout
from ...components.DropDown import DropDown

import scout


class PlotVarExplained(DashAction):
    def __init__(self, parent_cid: CID, dataset, input_pcy_cid: CID, select_plot_type_cid: CID):
        super().__init__(parent_cid, dataset)
        self.input_pcy_cid = input_pcy_cid
        self.select_plot_type_cid = select_plot_type_cid

    def plot(self, n_pcs, plot_type):
        var_explained = scout.ply.pca_explain_var(
            self.dataset.adata, layout=default_layout,
            plot_type=plot_type, n_pcs=n_pcs
        )
        return var_explained

    def setup_callbacks(self, app):
        output = Output(self.parent_cid.to_dict(), "figure")

        inputs = dict(
            n_pcs=Input(self.input_pcy_cid.to_dict(), "value"),
            plot_type=Input(self.select_plot_type_cid.to_dict(), "value"),
            _=Input("url", "pathname"),
        )

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(n_pcs, plot_type, _):
            return self.plot(n_pcs, plot_type)


class PlotCorrExplained(DashAction):
    def __init__(self, parent_cid: CID, dataset, input_pcy_cid: CID):
        super().__init__(parent_cid, dataset)
        self.input_pcy_cid = input_pcy_cid

    def plot(self, n_pcs):
        corr_explained = scout.ply.pca_explain_corr(
            self.dataset.adata, layout=default_layout, n_pcs=n_pcs
        )
        return corr_explained

    def setup_callbacks(self, app):
        output = Output(self.parent_cid.to_dict(), "figure")

        inputs = dict(
            n_pcs=Input(self.input_pcy_cid.to_dict(), "value"),
            _=Input("url", "pathname"),
        )
        @app.dash_app.callback(output=output, inputs=inputs)
        def _(n_pcs, _):
            return self.plot(n_pcs)



class CorrVarExplained(DashPlot):
    def __init__(self, dataset, page_id, loc_class):
        super().__init__(dataset, page_id, loc_class)
        n_pcs = self.dataset.adata.uns["pca"]["variance_ratio"].shape[0]
        self.children.update(
            input_pcy=InputField(
                cid=CID(self.page_id, self.loc_class, "input-n_pcs"), type="number",
                default=min(30, n_pcs), min=1, max=self.dataset.adata.uns["pca"]["variance_ratio"].shape[0],
            ),
            select_plot_type=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-plot_type"),
                options=["Bar", "Line", "Area", "Cumulative"], default="Cumulative"
            )
        )
        self.actions.update(
            plot_corr_explained=PlotCorrExplained(
                CID(self.page_id, self.loc_class, "subplot1"),
                self.dataset,
                self.children["input_pcy"].cid,
            ),
            plot_var_explained=PlotVarExplained(
                CID(self.page_id, self.loc_class, "subplot2"),
                self.dataset,
                self.children["input_pcy"].cid,
                self.children["select_plot_type"].cid,
            )
        )

    def create_layout(self) -> list:
        figure = html.Div([
            dcc.Loading(type="circle", children=[
                html.Div(
                    dcc.Graph(CID(self.page_id, self.loc_class, "subplot1").to_dict(), className=f"{self.loc_class}-left-plot")
                )
            ]),
            dcc.Loading(type="circle", children=[
                html.Div(
                    dcc.Graph(CID(self.page_id, self.loc_class, "subplot2").to_dict(), className=f"{self.loc_class}-right-plot")
                )
            ]),
        ], className=f"{self.loc_class}-body")

        return figure

    def get_sidebar_params(self) -> list:
        return [
            html.Div([
                html.Label("Plot Type"),
                self.children["select_plot_type"].create_layout(),
                self.children["select_plot_type"].get_stores()
            ], className="param-row-stacked"),
            html.Div([
                html.Label("Num. Components"),
                self.children["input_pcy"].create_layout(),
                self.children["input_pcy"].get_stores()
            ], className="param-row")
        ]
