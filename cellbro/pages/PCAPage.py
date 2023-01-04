import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

from ..components import components
from ..components.DashPage import DashPage
from ..util.DashAction import DashAction

from ..plots.PCA.CorrCircle import CorrCircle
from ..plots import projection as prj
from ..plots import PCA
from ..components.DashPlot import DashPlot
from ..components.Sidebar import Sidebar
from ..plots.PCA.CorrVarExplained import CorrVarExplained

import scout


class PCAPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.pca", "PCA", "pca", order)
        self.dataset = dataset
        self.components.update(
            pca_figure=prj.PCAProjection(self.dataset, self.page_id, "main"),
        )
        self.components.update(
            correlation_circle_figure=CorrCircle(
                dataset=self.dataset, page_id=self.page_id, loc_class="secondary",
                select_pcx_cid=self.components["pca_figure"].children["input_pcx"].cid,
                select_pcy_cid=self.components["pca_figure"].children["input_pcy"].cid
            ),
            bottom_figure=CorrVarExplained(self.dataset, self.page_id, "bottom")
        )

    def create_layout(self) -> list:
        self.components["left_sidebar"] = Sidebar(
            page_id=self.page_id, loc_class="main",
            title="PCA Projection Settings",
            params_children=self.components["pca_figure"].get_sidebar_params(),
            create_btn=False
        )


        self.components["bot_sidebar"] = Sidebar(
            page_id=self.page_id, loc_class="bottom",
            title="PCA Plots",
            params_children=self.components["bottom_figure"].get_sidebar_params(),
            create_btn=False
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
