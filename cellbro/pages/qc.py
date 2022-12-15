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

        inputs = {}
        for param in QC.qc_params.values():
            inputs[param.key] = Input(
                component_id=f"qc-{param.key}", component_property="value"
            )
        @app.dash_app.callback(output=output, inputs=inputs)
        def _(**kwargs):
            return self.apply(params=kwargs)

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


class SelectGeneListAction(DashAction):
    def apply(self, params):
        res = self.dataset.update_gene_lists(params["selected_gene"], params["gene_list"])
        return res, self.dataset.get_gene_lists()

    def setup_callbacks(self, app):
        @app.dash_app.callback(
            output=[
                Output("gene-list-dropdown", "value"),
                Output("gene-list-dropdown", "options"),
            ],
            inputs=[
                Input("gene-list-dropdown", "value"),
                Input("new-gene-list-button", "n_clicks"),
            ],
            state=[
                State("selected-gene", "children"),
                State("new-gene-list-input", "value"),
            ],
        )
        def _(gene_list, create_new_list, selected_gene, new_gene_list_name):
            if ctx.triggered_id == "new-gene-list-button":
                if create_new_list is None:
                    raise PreventUpdate
                if new_gene_list_name is None:
                    raise PreventUpdate
                if new_gene_list_name in self.dataset.get_gene_lists():
                    raise PreventUpdate
                self.dataset.adata.uns["gene_lists"][new_gene_list_name] = [selected_gene]
                return self.dataset.get_gene_lists(selected_gene), self.dataset.get_gene_lists()

            return self.apply(params=dict(selected_gene=selected_gene, gene_list=gene_list))



class QCPage(DashPage):
    def __init__(self, dataset, app, order):
        super().__init__("pages.qc", "QC", "/qc", order)
        self.dataset = dataset
        self.layout = self.create_layout()
        self.actions = dict(
            plot=PlotAction(dataset),
            filter=FilterAction(dataset),
            click=ClickAction(dataset),
            select_gene_list=SelectGeneListAction(dataset),
        )
        self.setup_callbacks(app)

    def create_layout(self):
        top_sidebar = html.Div(
            children=[
                html.Div(
                    [
                        html.H3("Filtering Settings"),
                    ],
                    className="sidebar-header",
                ),
                dcc.Loading(
                    type="circle",
                    children=[
                        html.Div(
                            children=[
                                QC.params_layout(),
                            ],
                            className="sidebar-parameters",
                        ),
                        html.Div(
                            [
                                dbc.Button(
                                    "Filter",
                                    color="primary",
                                    className="mr-1",
                                    id="filtering-submit",
                                ),
                            ],
                            className="sidebar-footer",
                        ),
                    ],
                ),
            ],
            className="top-sidebar sidebar",
        )

        bottom_sidebar = html.Div(
            children=[
                html.Div(
                    [
                        html.H3("QC Violin Plots"),
                    ],
                    className="sidebar-header",
                ),
                dcc.Loading(
                    type="circle",
                    children=[
                        # html.Div(children=[
                        # ], className="sidebar-parameters"),
                        # html.Div([
                        #     dbc.Button("Plot", color="primary", className="mr-1", id="qc-submit"),
                        # ], className="sidebar-footer")
                    ],
                ),
            ],
            id="qc-violin-sidebar",
            className="bottom-sidebar sidebar",
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
            html.Div(
                id="top",
                className="top",
                children=[top_sidebar, main_figure, secondary_figure],
            ),
            html.Div(
                id="bottom", className="bottom", children=[bottom_sidebar, bottom_figure]
            ),
        ]
        return layout