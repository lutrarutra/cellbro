import scanpy as sc
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, State, ctx, dcc, html
from dash.exceptions import PreventUpdate

import cellbro.plots.QC as QC
from cellbro.util.DashPage import DashPage
from cellbro.util.DashAction import DashAction
import cellbro.util.Components as Components

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
        for param in QC.qc_params.values():
            output.append(
                Output(component_id=f"{self.page_id_prefix}-{param.key}", component_property="value")
            )

        inputs = {
            "submit": Input(
                component_id=f"{self.page_id_prefix}-filtering-submit", component_property="n_clicks"
            ),
        }
        state = {}
        for param in QC.qc_params.values():
            state[param.key] = State(
                component_id=f"{self.page_id_prefix}-{param.key}", component_property="value"
            )

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(**kwargs):
            return self.apply(params=kwargs)


class ClickAction(DashAction):
    def apply(self, params):
        data = params["clickData"]["points"][0]
        gene = data["hovertext"]
        element = Components.create_gene_card(gene, self.dataset)
        return [element]

    def setup_callbacks(self, app):
        @app.dash_app.callback(
            [Output(f"{self.page_id_prefix}-dispersion-info", "children")], [Input(f"{self.page_id_prefix}-dispersion-plot", "clickData")]
        )
        def _(clickData):
            if clickData is None:
                raise PreventUpdate
            # return QC.on_hover(clickData["points"][0], self.dataset)
            return self.apply(params=dict(clickData=clickData))


class PerformQC(DashAction):
    def apply(self):
        if not "pct_counts_mt" in self.dataset.adata.obs.columns:
            QC.apply_mt_qc(self.dataset)

        if not "cv2" in self.dataset.adata.var.columns or not "mu" in self.dataset.adata.var.columns:
            QC.apply_dispersion_qc(self.dataset)

        return [{"display": "none"}, {"display": "block"}]

    def plot(self, params):
        return QC.plot(self.dataset, params)

    def setup_callbacks(self, app):
        state = {}
        for param in QC.qc_params.values():
            state[param.key] = State(
                component_id=f"{self.page_id_prefix}-{param.key}", component_property="value"
            )

        @app.dash_app.callback(
            output=[
                Output(f"{self.page_id_prefix}-top-sidebar-temp", "style"),
                Output(f"{self.page_id_prefix}-top-sidebar", "style"),
                Output(component_id=f"{self.page_id_prefix}-mt-plot", component_property="figure"),
                Output(component_id=f"{self.page_id_prefix}-dispersion-plot", component_property="figure"),
                Output(component_id=f"{self.page_id_prefix}-violin-plot", component_property="figure"),
            ],
            inputs=dict(submit=Input(f"{self.page_id_prefix}-apply-btn", "n_clicks")),
            state=state,
        )
        def _(submit, **kwargs):
            if submit:
                return self.apply() + self.plot(params=kwargs)

            if self.dataset.qc_done():
                return [{"display": "none"}, {"display": "block"}] + self.plot(params=kwargs)

            raise PreventUpdate
            # return [{"display": "block"}, {"display": "none"}, None, None, None]

class QCPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.qc", "QC", "qc", order)
        self.dataset = dataset
        self.actions.update(
            filter=FilterAction(dataset, self.id),
            click=ClickAction(dataset, self.id),
            top_sidebar_temp=Components.HideSidebar(page_id_prefix=self.id, id=f"{self.id}-top-sidebar-temp"),
            perform_qc=PerformQC(dataset, self.id),
        )

    def create_layout(self):

        top_sidebar = Components.create_sidebar(
            id=f"{self.id}-top-sidebar", class_name="top-sidebar",
            title="Quality Control Parameters",
            params_children=self._filter_params_layout(),
            btn_id=f"{self.id}-filtering-submit", btn_text="Filter"
        )

        temp_top_sidebar = Components.create_sidebar(
            id=f"{self.id}-top-sidebar-temp", class_name="top-sidebar",
            title="Perform Quality Control",
            params_children=self._qc_params_layout(),
            btn_id=f"{self.id}-apply-btn", btn_text="Quality Control"
        )

        bot_sidebar = Components.create_sidebar(
            id=f"{self.id}-bot-sidebar", class_name="bot-sidebar",
            title="Quality Control Violin Plots",
            params_children=[],
            btn_id=None, btn_text="Filter"
        )

        main_figure = html.Div(
            children=[
                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(dcc.Graph(id=f"{self.id}-mt-plot", className="main-plot"))
                            ],
                        )
                    ],
                    id=f"{self.id}-mt-figure",
                    className="main-figure",
                )
            ],
            className="main",
        )

        secondary_figure = html.Div(
            children=[
                html.Div(
                    id=f"{self.id}-dispersion-info",
                    children=[
                        Components.create_gene_card(None, self.dataset),
                    ],
                    className="secondary-select top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id=f"{self.id}-dispersion-plot", className="secondary-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    id=f"{self.id}-dispersion-figure",
                    className="secondary-figure",
                ),
            ],
            className="secondary",
        )

        bottom_figure = html.Div(
            children=[
                dcc.Loading(
                    type="circle",
                    children=[
                        html.Div(
                            dcc.Graph(id=f"{self.id}-violin-plot", className="bottom-plot")
                        )
                    ],
                )
            ],
            id=f"{self.id}-violin-figure",
            className="bottom-figure",
        )

        layout = [
            html.Div(
                className="top",
                children=[top_sidebar, temp_top_sidebar, main_figure, secondary_figure],
            ),
            html.Div(
                className="bottom", children=[bot_sidebar, bottom_figure]
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
                    id=f"{self.id}-{key}", type=param.input_type, value=param.value,
                    step=param.step if param.step != None else 0.1, className="param-input",
                ),
            ], className="param-row") for key, param in QC.qc_params.items()
        ]

    def _qc_params_layout(self):
        return [
            html.Div(children=[
                html.Label(
                    "Perform QC before filtering", className="param-label",
                ),
            ], className="param-row")
        ]


