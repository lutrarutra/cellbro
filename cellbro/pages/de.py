import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

import scout

import cellbro.plots.DE as DE
import cellbro.util.Components as Components
from cellbro.util.DashAction import DashAction
from cellbro.util.DashPage import DashPage


class ApplyDE(DashAction):
    def apply(self, params):
        return DE.apply(self.dataset, params=params) + [{"new": True}]

    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-volcano-groupby", component_property="options"),
            Output(component_id=f"{self.page_id_prefix}-volcano-groupby", component_property="value"),
            Output(component_id=f"{self.page_id_prefix}-volcano-reference", component_property="options"),
            Output(component_id=f"{self.page_id_prefix}-volcano-reference", component_property="value"),
            Output(
                component_id=f"{self.page_id_prefix}-volcano-groupby-container", component_property="style"
            ),
            Output(
                component_id=f"{self.page_id_prefix}-volcano-reference-container", component_property="style"
            ),
            Output(component_id="groupby-store", component_property="data"),
        ]
        inputs = {
            "submit": Input(component_id=f"{self.page_id_prefix}-submit", component_property="n_clicks"),
        }
        state = {
            "groupby": State(component_id=f"{self.page_id_prefix}-groupby", component_property="value"),
        }
        for param in DE.de_params.values():
            state[param.key] = State(
                component_id=f"de-{param.key}", component_property="value"
            )

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(**kwargs):
            return self.apply(params=kwargs)


class PlotVolcano(DashAction):
    def apply(self, params):
        return [DE.plot_de_volcano(dataset=self.dataset, params=params)]

    def setup_callbacks(self, app):
        outputs = [
            Output(component_id=f"{self.page_id_prefix}-volcano-plot", component_property="figure"),
        ]
        inputs = {
            "groupby": Input(
                component_id=f"{self.page_id_prefix}-volcano-groupby", component_property="value"
            ),
            "reference": Input(
                component_id=f"{self.page_id_prefix}-volcano-reference", component_property="value"
            ),
        }

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(**kwargs):
            if kwargs["groupby"] is None:
                raise PreventUpdate
            return self.apply(params=kwargs)


class PlotPvalHistogram(DashAction):
    def apply(self, params):
        return [DE.plot_pval_histogram(dataset=self.dataset, params=params)]

    def setup_callbacks(self, app):
        outputs = [
            Output(component_id=f"{self.page_id_prefix}-de-secondary-plot", component_property="figure"),
        ]
        inputs = {
            "groupby": Input(
                component_id=f"{self.page_id_prefix}-volcano-groupby", component_property="value"
            ),
            "reference": Input(
                component_id=f"{self.page_id_prefix}-volcano-reference", component_property="value"
            ),
        }

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(**kwargs):
            if kwargs["groupby"] is None:
                raise PreventUpdate
            return self.apply(params=kwargs)


class ClickAction(DashAction):
    def apply(self, params):
        gene = params["click_data"]["points"][0]["hovertext"]
        element = Components.create_gene_card(gene, self.dataset)
        return [element]

    def setup_callbacks(self, app):
        outputs = [Output(f"{self.page_id_prefix}-volcano-info", "children")]
        inputs = {
            "click_data": Input(f"{self.page_id_prefix}-volcano-plot", "clickData"),
        }

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(**kwargs):
            if kwargs["click_data"] is None:
                raise PreventUpdate
            return self.apply(params=kwargs)


class DEPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.de", "DE", "de", order)
        self.dataset = dataset
        self.actions.update(
            apply_de=ApplyDE(dataset=self.dataset, page_id_prefix=self.id),
            plot_volcano=PlotVolcano(dataset=self.dataset, page_id_prefix=self.id),
            plot_pval_histogram=PlotPvalHistogram(dataset=self.dataset, page_id_prefix=self.id),
            click_action=ClickAction(dataset=self.dataset, page_id_prefix=self.id),
        )

    def create_layout(self):
        top_sidebar = Components.create_sidebar(
            id=f"{self.id}-top-sidebar", btn_id=f"{self.id}-submit",
            title="Differential Expression Settings",
            params_children=self._params_layout(),
            class_name="top-sidebar"
        )

        main = html.Div(
            children=[
                html.Div(
                    children=[
                        Components.create_gene_card(None, self.dataset),
                    ],
                    id=f"{self.id}-volcano-info",
                    className="main-select top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id=f"{self.id}-volcano-plot", className="main-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    id=f"{self.id}-volcano-figure",
                    className="main-figure",
                ),
            ],
            className="main",
        )

        secondary_options = {"pval_histogram": "P-Value Histogram"}

        secondary = html.Div(
            children=[
                html.Div(
                    children=[
                        # Volcano Group By Select
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Label("Plot Type"),
                                dcc.Dropdown(
                                    options=secondary_options,
                                    value="pval_histogram",
                                    id=f"{self.id}-secondary-type",
                                    clearable=False,
                                ),
                            ],
                            className="param-column",
                        ),
                    ],
                    id=f"{self.id}-de-secondary-select",
                    className="main-select top-parameters",
                ),
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
                    id=f"{self.id}-de-secondary-figure",
                    className="secondary-figure",
                ),
            ],
            className="secondary",
        )

        bot_sidebar = Components.create_sidebar(
            id=f"{self.id}-bot-sidebar", class_name="bot-sidebar",
            title="Empty", 
            params_children=[],
        )

        layout = [
            html.Div(
                className="top", children=[top_sidebar, main, secondary]
            ),
            html.Div(
                className="bottom",
                children=[
                    # bottom_sidebar, bottom_figure
                    bot_sidebar
                ],
            ),
        ]
        return layout

    def _params_layout(self):
        cats = self.dataset.get_categoric()
        divs = []

        groups = DE.get_groupby_options(dataset=self.dataset)
        refs = DE.get_reference_options(
            dataset=self.dataset, groupby=next(iter(groups), None)
        )
        # Volcano Group By Select
        divs.append(
            html.Div(
                id=f"{self.id}-volcano-groupby-container",
                children=[
                    html.Label("Volcano GroupBy"),
                    html.Div(
                        [
                            dcc.Dropdown(
                                options=groups,
                                value=next(iter(groups), None),
                                id=f"{self.id}-volcano-groupby",
                                clearable=False,
                            ),
                        ],
                        className="param-select",
                    ),
                ],
                className="param-row-stacked",
                style={"display": "none" if len(groups) == 0 else "block"},
            )
        )

        # Volcano Reference Select i.e. 'KO vs. Rest'
        divs.append(
            html.Div(
                id=f"{self.id}-volcano-reference-container",
                children=[
                    html.Label("Volcano Reference"),
                    html.Div(
                        [
                            dcc.Dropdown(
                                options=refs,
                                value=next(iter(refs), None),
                                id=f"{self.id}-volcano-reference",
                                clearable=False,
                            ),
                        ],
                        className="param-select",
                    ),
                ],
                className="param-row-stacked",
                style={"display": "none" if len(groups) == 0 else "block"},
            )
        )

        divs.append(
            html.Div(
                children=[
                    html.Label(
                        "Group By",
                        className="param-label",
                    ),
                    html.Div(
                        [
                            dcc.Dropdown(
                                options=cats,
                                value=cats[0],
                                id=f"de-groupby",
                                clearable=False,
                            )
                        ],
                        className="param-select",
                    ),
                ],
                className="param-row-stacked",
            )
        )

        for key, param in DE.de_params.items():
            divs.append(
                html.Div(
                    children=[
                        html.Label(
                            param.name,
                            className="param-label",
                        ),
                        html.Div(
                            [
                                dcc.Dropdown(
                                    id=f"de-{key}",
                                    value=param.default,
                                    options=param.allowed_values,
                                    clearable=False,
                                )
                            ],
                            className="param-select",
                        ),
                    ],
                    className="param-row-stacked",
                )
            )

        # layout = html.Div(children=divs)
        return divs
