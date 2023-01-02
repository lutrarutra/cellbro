import scanpy as sc
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, State, ctx, dcc, html
from dash.exceptions import PreventUpdate

from ..components.DashPage import DashPage
from ..util.DashAction import DashAction
from ..components.Sidebar import Sidebar
from ..plots import QC
from ..components.CID import CID, LocClass

class FilterAction(DashAction):
    def apply(self, params):
        # Makes sure that filtering is not done on initial load
        submit = params.pop("submit")
        if submit is None:
            return list(params.values())

        self.dataset.adata = self.dataset.adata[
            self.dataset.adata.obs.pct_counts_mt < self.params["pct_counts_mt"], :
        ].copy()
        sc.pp.filter_cells(self.dataset.adata, min_genes=self.params["min_genes"])
        sc.pp.filter_genes(self.dataset.adata, min_cells=self.params["min_cells"])

        return list(self.params.values())

    def setup_callbacks(self, app):
        output = []
        for param in QC.qc_tools.qc_params.values():
            output.append(
                Output(f"{self.page_id}-{param.key}", "value")
            )

        inputs = dict(
            submit=Input(f"{self.page_id}-main-sidebar-apply_btn", "n_clicks"),
        )

        for param in QC.qc_tools.qc_params.values():
            inputs[param.key] = State(f"{self.page_id}-{param.key}", "value")

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(**kwargs):
            return self.apply(params=kwargs)


class QCPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.qc", "QC", "qc", order)
        self.dataset = dataset
        self.actions.update(
            filter=FilterAction(CID(self.page_id, LocClass.static, "filter_action"), dataset),
        )
        self.components.update(
            main_figure=QC.MTPlot(self.dataset, self.page_id, "main"),
            secondary_figure=QC.DispersionPlot(self.dataset, self.page_id, "secondary"),
            bot_figure=QC.QCViolins(self.dataset, self.page_id, "bottom"),
        )

    def create_layout(self):

        self.components["left_sidebar"] = Sidebar(
            page_id=self.page_id, loc_class="main",
            title="Quality Control Parameters",
            params_children=self._filter_params_layout(),
            create_btn=True, btn_text="Filter"
        )

        # Something smarter
        self.components["right_sidebar"] = Sidebar(
            page_id=self.page_id, loc_class="secondary",
            title="Quality Control Parameters",
            params_children=self._qc_params_layout(),
            create_btn=True, btn_text="Filter"
        )

        self.components["bot_sidebar"] = Sidebar(
            page_id=self.page_id, loc_class="bottom",
            title="Quality Control Violin Plots",
            params_children=[],
        )

        layout = [
            html.Div(
                className="top",
                children=[
                    self.components["left_sidebar"].create_layout(),
                    self.components["main_figure"].create_layout(),
                    self.components["secondary_figure"].create_layout(),
                    self.components["right_sidebar"].create_layout()
                ],
            ),
            html.Div(
                className="bottom", children=[
                    self.components["bot_sidebar"].create_layout(),
                    self.components["bot_figure"].create_layout(),
                ]
            ),
        ]
        return layout

    def _filter_params_layout(self):
        return [
            html.Div(children=[
                html.Label(
                    param.name, className="param-label",
                ),
                dcc.Input(
                    id=f"{self.page_id}-{key}", type=param.input_type, value=param.value,
                    step=param.step if param.step != None else 0.1, className="param-input",
                ),
            ], className="param-row") for key, param in QC.qc_tools.qc_params.items()
        ]

    def _qc_params_layout(self):
        return [
            html.Div(children=[
                html.Label(
                    "Perform QC before filtering", className="param-label",
                ),
            ], className="param-row")
        ]


