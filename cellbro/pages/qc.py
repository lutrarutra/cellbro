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


# class DispersionPlot(DashAction):
#     def apply(self, params):
#         fig = px.scatter(
#             self.dataset.adata.var.reset_index(),
#             x="mu",
#             y="cv2",
#             log_x=True,
#             log_y=True,
#             color_continuous_scale=px.colors.sequential.Viridis,
#             hover_name="index",
#         )
#         fig.update_traces(
#             marker=dict(size=5, line=dict(width=1, color="DarkSlateGrey"))
#         )
#         fig.update_layout(QC.figure_layout)
#         fig.update_layout(xaxis_title="Log Mean Expression", yaxis_title="CV^2")
#         return fig

#     def setup_callbacks(self, dash_app):
#         output = [Output(component_id="dispersion-plot", component_property="figure")]
#         inputs = {}
#         for param in QC.qc_params.values():
#             inputs[param.key] = Input(
#                 component_id=f"qc-{param.key}", component_property="value"
#         )

class PlotAction(DashAction):
    def apply(self, params):
        return QC.plot(self.dataset, params)

    def setup_callbacks(self, app):
        output = [
            Output(component_id="mt-plot", component_property="figure"),
            Output(component_id="dispersion-plot", component_property="figure"),
            Output(component_id="qc-violin-plot", component_property="figure"),
        ]

        inputs = {
            "qc_store": Input("qc-store", "data")
        }
        for param in QC.qc_params.values():
            inputs[param.key] = Input(
                component_id=f"qc-{param.key}", component_property="value"
            )
        @app.dash_app.callback(output=output, inputs=inputs)
        def _(qc_store, **kwargs):
            if self.dataset.qc_done() or qc_store["qc_done"]:
                return self.apply(params=kwargs)

            raise PreventUpdate

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
                Output(component_id=f"qc-{param.key}", component_property="value")
            )

        inputs = {
            "submit": Input(
                component_id="filtering-submit", component_property="n_clicks"
            ),
        }
        state = {}
        for param in QC.qc_params.values():
            state[param.key] = State(
                component_id=f"qc-{param.key}", component_property="value"
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
            [Output("dispersion-info", "children")], [Input("dispersion-plot", "clickData")]
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

        return [{"display": "none"}, {"display": "block"}, {"qc_done": True}]

    def setup_callbacks(self, app):
        @app.dash_app.callback(
            output=[
                Output("qc-top-sidebar-temp", "style"),
                Output("qc-top-sidebar", "style"),
                Output("qc-store", "data")
            ],
            inputs=[
                Input("perform-qc-btn", "n_clicks"),
            ],
        )
        def _(submit):
            if submit:
                return self.apply()

            if self.dataset.qc_done():
                return [{"display": "none"}, {"display": "block"}, {"qc_done": True}]

            return [{"display": "block"}, {"display": "none"}, {"qc_done": False}]

class QCPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.qc", "QC", "/qc", order)
        self.dataset = dataset
        self.actions.update(
            plot=PlotAction(dataset),
            filter=FilterAction(dataset),
            click=ClickAction(dataset),
            top_sidebar_temp=Components.HideSidebar(id=f"qc-top-sidebar-temp"),
            perform_qc=PerformQC(dataset),
        )

    def create_layout(self):

        top_sidebar = Components.create_sidebar(
            id="qc-top-sidebar", class_name="top-sidebar",
            title="Quality Control Parameters",
            params_children=self._filter_params_layout(),
            btn_id="filtering-submit", btn_text="Filter"
        )

        temp_top_sidebar = Components.create_sidebar(
            id="qc-top-sidebar-temp", class_name="top-sidebar",
            title="Perform Quality Control",
            params_children=self._qc_params_layout(),
            btn_id="perform-qc-btn", btn_text="Quality Control"
        )

        bot_sidebar = Components.create_sidebar(
            id="qc-bot-sidebar", class_name="bot-sidebar",
            title="Quality Control Violin Plots",
            params_children=[],
            btn_id=None, btn_text="Filter"
        )

        main_figure = html.Div(
            children=[
                html.Div(
                    [
                        dcc.Loading(
                            id="loading-mt",
                            type="circle",
                            children=[
                                html.Div(dcc.Graph(id="mt-plot", className="main-plot"))
                            ],
                        )
                    ],
                    id="mt-figure",
                    className="main-figure",
                )
            ],
            className="main",
        )

        secondary_figure = html.Div(
            children=[
                html.Div(
                    id="dispersion-info",
                    children=[
                        Components.create_gene_card(None, self.dataset),
                    ],
                    className="secondary-select top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            id="loading-dispersion",
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id="dispersion-plot", className="secondary-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    id="dispersion-figure",
                    className="secondary-figure",
                ),
            ],
            className="secondary",
        )

        bottom_figure = html.Div(
            children=[
                dcc.Loading(
                    id="loading-qc-violin",
                    className="loading-bottom",
                    type="circle",
                    children=[
                        html.Div(
                            dcc.Graph(id="qc-violin-plot", className="bottom-plot")
                        )
                    ],
                )
            ],
            id="qc-violin-figure",
            className="bottom-figure",
        )

        layout = [
            dcc.Store(id="qc-store"),
            html.Div(
                id="top",
                className="top",
                children=[top_sidebar, temp_top_sidebar, main_figure, secondary_figure],
            ),
            html.Div(
                id="bottom", className="bottom", children=[bot_sidebar, bottom_figure]
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
                    id=f"qc-{key}", type=param.input_type, value=param.value,
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


