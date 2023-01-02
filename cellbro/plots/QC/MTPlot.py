from dash import Input, Output, State, ctx, dcc, html
from dash.exceptions import PreventUpdate

from ...components.DashPlot import DashPlot
from ...util.DashAction import DashAction
from .qc_tools import default_layout

import scout


class Plot(DashAction):
    def plot(self):
        fig = scout.ply.mt_plot(
            self.dataset.adata, layout=default_layout
        )
        return fig

    def setup_callbacks(self, app):
        output = Output(self.parent_cid.to_str(), "figure")
        inputs = dict(
            submit=Input(f"{self.page_id}-{self.loc_class}-sidebar-apply_btn", "n_clicks"),
        )

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(submit):
            if self.dataset.qc_done():
                return self.plot()

            raise PreventUpdate


class MTPlot(DashPlot):
    def __init__(self, dataset, page_id, loc_class):
        super().__init__(dataset, page_id, loc_class)
        self.actions.update(
            plot=Plot(self.cid, dataset)
        )

    def create_layout(self):
        figure = html.Div([
            html.Div([
                dcc.Loading(type="circle", children=[
                    html.Div(
                        dcc.Graph(id=f"{self.page_id}-{self.loc_class}-plot", className=f"{self.loc_class}-plot")
                    )
                ])
            ], className=f"{self.loc_class}-body")
        ], className=f"{self.loc_class}")

        return figure
