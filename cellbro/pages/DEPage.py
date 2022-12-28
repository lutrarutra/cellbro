import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

import scout

from ..components import components
from ..util.DashAction import DashAction
from ..components.DashPage import DashPage

from ..plots import DE

class PlotPvalHistogram(DashAction):
    def apply(self, params):
        return [DE.de_tools.plot_pval_histogram(dataset=self.dataset, params=params)]

    def setup_callbacks(self, app):
        outputs = [
            Output(component_id=f"{self.page_id_prefix}-de-secondary-plot", component_property="figure"),
        ]
        inputs = {
            "groupby": Input(
                component_id=f"{self.page_id_prefix}-main-groupby", component_property="value"
            ),
            "reference": Input(
                component_id=f"{self.page_id_prefix}-main-reference", component_property="value"
            ),
        }

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(**kwargs):
            if kwargs["groupby"] is None:
                raise PreventUpdate
            return self.apply(params=kwargs)

class ApplyDE(DashAction):
    def apply(self, groupby, reference, params):
        scout.tl.rank_marker_genes(self.dataset.adata, groupby=groupby, reference=reference, **params)

        return dict(update=True, groupby=groupby)

    def setup_callbacks(self, app):
        output = Output("de-store", "data")

        inputs = dict(submit=Input(component_id=f"{self.page_id_prefix}-submit", component_property="n_clicks"))

        state = dict(
            groupby=State(component_id=f"{self.page_id_prefix}-groupby", component_property="value"),
            reference=State(component_id=f"{self.page_id_prefix}-reference", component_property="value")
        )

        for param in DE.de_tools.de_params.values():
            state[param.key] = State(
                component_id=f"de-{param.key}", component_property="value"
            )

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, groupby, reference, **kwargs):
            if submit is None:
                raise PreventUpdate

            return self.apply(groupby, reference, params=kwargs)

        # Update References
        output = Output(f"{self.page_id_prefix}-reference", "options")
        inputs = dict(
            groupby=Input(f"{self.page_id_prefix}-groupby", "value"),
        )

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(groupby):
            cats = self.dataset.get_categoric()

            if len(cats) == 0:
                ref_options = []
            else:
                ref_options = self.dataset.adata.obs[groupby].cat.categories.to_list()

            return ["rest"] + ref_options

    def get_sidebar_params(self) -> list:
        cats = self.dataset.get_categoric()

        divs = []

        divs.append(
            html.Div(
                children=[
                    html.Label("Group By", className="param-label",),
                    html.Div([
                        dcc.Dropdown(options=cats, value=cats[0], id=f"{self.page_id_prefix}-groupby", clearable=False)
                    ], className="param-select"),
                ],
                className="param-row-stacked",
            )
        )

        divs.append(
            html.Div(
                children=[
                    html.Label("Reference", className="param-label"),
                    html.Div([
                        dcc.Dropdown(
                            options=["rest"], value="rest", id=f"{self.page_id_prefix}-reference", clearable=False
                        )
                    ], className="param-select"),
                ],
                className="param-row-stacked",
            )
        )

        for key, param in DE.de_tools.de_params.items():
            divs.append(
                html.Div(
                    children=[
                        html.Label(param.name, className="param-label",),
                        html.Div([
                            dcc.Dropdown(
                                id=f"{self.page_id_prefix}-{key}", value=param.default, options=param.allowed_values,
                                clearable=False,
                            )
                        ], className="param-select"),
                    ],
                    className="param-row-stacked",
                )
            )

        # layout = html.Div(children=divs)
        return divs


class DEPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.de", "DE", "de", order)
        self.dataset = dataset
        self.actions.update(
            de_apply=ApplyDE(dataset=self.dataset, page_id_prefix=self.id),
            plot_pval_histogram=PlotPvalHistogram(dataset=self.dataset, page_id_prefix=self.id),
        )
        self.components.update(
            main_figure=DE.DEVolcano(dataset=self.dataset, page_id_prefix=self.id, loc_class="main"),
        )

    def create_layout(self):
        self.components["top_sidebar"] = components.Sidebar(
            page_id_prefix=self.id, apply_btn_id=f"{self.id}-submit",
            title="Differential Expression Settings",
            params_children=self.actions["de_apply"].get_sidebar_params(),
            row="top", side="left",
        )

        plot_type_params = components.FigureHeaderTab(self.id, tab_label="Type", children=[
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

        figure_params = components.FigureHeader(self.id, tabs=[plot_type_params])

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

        self.components["bot_sidebar"] = components.Sidebar(
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
