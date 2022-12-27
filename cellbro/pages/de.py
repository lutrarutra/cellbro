import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

import scout

import cellbro.plots.DE as DE
import cellbro.util.Components as Components
from cellbro.util.DashAction import DashAction
from cellbro.util.DashPage import DashPage
from ..plots.DEVolcanoFig import DEVolcanoFig

class PlotPvalHistogram(DashAction):
    def apply(self, params):
        return [DE.plot_pval_histogram(dataset=self.dataset, params=params)]

    def setup_callbacks(self, app):
        outputs = [
            Output(component_id=f"{self.page_id_prefix}-de-secondary-plot", component_property="figure"),
        ]
        inputs = {
            "groupby": Input(
                component_id=f"{self.page_id_prefix}-main-volcano-groupby", component_property="value"
            ),
            "reference": Input(
                component_id=f"{self.page_id_prefix}-main-volcano-reference", component_property="value"
            ),
        }

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(**kwargs):
            if kwargs["groupby"] is None:
                raise PreventUpdate
            return self.apply(params=kwargs)


class ApplyDE(DashAction):
    def apply(self, params):
        groupby = params["groupby"]

        scout.tl.rank_marker_genes(self.dataset.adata, groupby=groupby)

        return dict(update=True, groupby=groupby)

    def setup_callbacks(self, app):
        output = Output("de-store", "data")

        inputs = dict(submit=Input(component_id=f"{self.page_id_prefix}-submit", component_property="n_clicks"))

        state = {
            "groupby": State(component_id=f"{self.page_id_prefix}-groupby", component_property="value"),
        }

        for param in DE.de_params.values():
            state[param.key] = State(
                component_id=f"de-{param.key}", component_property="value"
            )

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, **kwargs):
            if submit is None:
                raise PreventUpdate

            return self.apply(params=kwargs)


class DEPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.de", "DE", "de", order)
        self.dataset = dataset
        self.actions.update(
            de_apply=ApplyDE(dataset=self.dataset, page_id_prefix=self.id),
            plot_pval_histogram=PlotPvalHistogram(dataset=self.dataset, page_id_prefix=self.id),
        )
        self.components.update(
            main_figure=DEVolcanoFig(dataset=self.dataset, page_id_prefix=self.id, loc_class="main"),
        )

    def create_layout(self):
        self.components["top_sidebar"] = Components.Sidebar(
            page_id_prefix=self.id, apply_btn_id=f"{self.id}-submit",
            title="Differential Expression Settings",
            params_children=self.components["main_figure"].get_sidebar_params(),
            row="top", side="left",
        )

        plot_type_params = Components.FigureParamTab(self.id, tab_label="Type", children=[
            # Volcano Group By Select
            html.Div([
                html.Label("Plot Type"),
                dcc.Dropdown(
                    options={"pval_histogram": "P-Value Histogram"},
                    value="pval_histogram",
                    id=f"{self.id}-secondary-type",
                    clearable=False,
                ),
            ], className="param-row-stacked")
        ])

        figure_params = Components.FigureParams(self.id, tabs=[plot_type_params])

        secondary = html.Div(
            children=[
                html.Div(figure_params.create_layout(), className="fig-header"),
                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id=f"{self.id}-de-secondary-plot",
                                        className="secondary-plot",
                                    )
                                )
                            ],
                        )
                    ],
                    id=f"{self.id}-de-secondary-body",
                    className="secondary-body",
                ),
            ],
            className="secondary",
        )

        self.components["bot_sidebar"] = Components.Sidebar(
            page_id_prefix=self.id,  row="bot", side="left",
            title="Empty",
            params_children=[], apply_btn_id=None
        )

        layout = [
            html.Div(
                className="top", children=[
                    self.components["top_sidebar"].create_layout(), 
                    self.components["main_figure"].create_layout(), secondary
                ]
            ),
            html.Div(
                className="bottom",
                children=[
                    self.components["bot_sidebar"].create_layout()
                ],
            ),
        ]
        return layout
