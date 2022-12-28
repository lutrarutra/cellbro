import dash
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, State, dcc, html, ctx
from dash.exceptions import PreventUpdate

from ..components import components
from ..components.DashPage import DashPage
from ..plots.Heatmap import Heatmap
from ..plots.Violin import Violin
from ..plots import projection as prj

import scout

class CellsPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.cells", "Cells", "cells", order)
        self.dataset = dataset
        self.components.update(
            projection=prj.Projection(self.dataset, self.id, "main"),
            violin=Violin(self.dataset, self.id, "secondary"),
            heatmap=Heatmap(self.dataset, self.id, "bottom"),
        )

    def create_layout(self) -> list:
        self.components["left_sidebar"] = components.Sidebar(
            page_id_prefix=self.id, row="top", side="left",
            title="Projection Settings",
            params_children=self.components["projection"].get_sidebar_params(),
            apply_btn_id=f"{self.id}-projection-submit", btn_text="Apply Projection"
        )

        main_figure = self.components["projection"].create_layout()

        self.components["bot_sidebar"] = components.Sidebar(
            page_id_prefix=self.id, row="bot", side="left",
            title="Heatmap Settings",
            params_children=self.components["heatmap"].get_sidebar_params(),
            apply_btn_id=f"{self.id}-heatmap-submit", btn_text="Plot"
        )

        bottom_figure = self.components["heatmap"].create_layout()
        violin_layout = self.components["violin"].create_layout()

        return [
            html.Div(
                className="top",
                children=[self.components["left_sidebar"].create_layout(), main_figure, violin_layout],
            ),
            html.Div(
                className="bottom",
                children=[self.components["bot_sidebar"].create_layout(), bottom_figure],
            ),
        ]
