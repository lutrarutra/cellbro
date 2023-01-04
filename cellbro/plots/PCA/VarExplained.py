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

import scout


class PlotVarExplained(DashAction):
    def __init__(self, parent_cid: CID, dataset, input_pcy_cid: CID):
        super().__init__(parent_cid, dataset)
        self.input_pcy_cid = input_pcy_cid

    def plot(self, plot_type, n_pcs):
        fig = scout.ply.pca_explain_var(
            self.dataset.adata, layout=default_layout,
            plot_type=plot_type, n_pcs=n_pcs
        )

        return fig

    def setup_callbacks(self, app):
        output = Output(self.page_id.to_dict(), "figure")

        inputs = dict(
            n_pcs=Input(self.input_pcy_cid.to_dict(), "value"),
            plot_type=Input(f"{self.page_id}-hist-plot_type", "value"),
        )

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(plot_type, n_pcs):
            return self.plot(plot_type, n_pcs)
