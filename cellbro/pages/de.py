import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

import cellbro.plots.DE as DE
import scout
from cellbro.plots.DashPage import DashPage
from cellbro.util.Components import create_gene_card


class DEPage(DashPage):
    def __init__(self, dataset, dash_app):
        super().__init__("pages.de", "DE", "/de", 4)
        self.dataset = dataset
        self.layout = self.create_layout()
        self.setup_callbacks(dash_app)

    def plot(self, **kwargs) -> list:
        volcano = DE.plot_de_volcano(dataset=self.dataset, params=kwargs)
        secondary = DE.plot_pval_histogram(dataset=self.dataset, params=kwargs)
        return [volcano, secondary]

    def apply(self, **kwargs):
        return DE.apply(dataset=self.dataset, params=kwargs)

    def setup_callbacks(self, dash_app):
        @dash_app.callback(**self._get_apply_callbacks())
        def _(**kwargs):
            return self.apply(**kwargs)

        @dash_app.callback(**self._get_plot_callbacks())
        def _(**kwargs):
            if kwargs["groupby"] is None:
                raise PreventUpdate
            return self.plot(**kwargs)

        @dash_app.callback(**self._get_on_click_callbacks())
        def _(click_data):
            if click_data is None:
                raise PreventUpdate
            gene = click_data["points"][0]["hovertext"]
            element = create_gene_card(gene, self.dataset)
            return [element]

    def create_layout(self):
        top_sidebar = html.Div(
            children=[
                html.Div(
                    [
                        html.H3("Differential Expression Settings"),
                    ],
                    className="sidebar-header",
                ),
                dcc.Loading(
                    type="circle",
                    children=[
                        html.Div(
                            children=[
                                self._params_layout(),
                            ],
                            className="sidebar-parameters",
                        ),
                        html.Div(
                            [
                                dbc.Button(
                                    "Apply",
                                    color="primary",
                                    className="mr-1",
                                    id="de-submit",
                                ),
                            ],
                            className="sidebar-footer",
                        ),
                    ],
                ),
            ],
            className="top-sidebar sidebar",
        )

        main = html.Div(
            children=[
                html.Div(
                    children=[
                        create_gene_card(None, self.dataset),
                    ],
                    id="de-volcano-info",
                    className="main-select top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            id="loading-de-volcano",
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id="de-volcano-plot", className="main-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    id="de-volcano-figure",
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
                                    id="secondary-type",
                                    clearable=False,
                                ),
                            ],
                            className="param-column",
                        ),
                    ],
                    id="de-secondary-select",
                    className="main-select top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            id="loading-de-secondary",
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id="de-secondary-plot",
                                        className="secondary-plot",
                                    )
                                )
                            ],
                        )
                    ],
                    id="de-secondary-figure",
                    className="secondary-figure",
                ),
            ],
            className="secondary",
        )

        layout = [
            html.Div(
                id="top", className="top", children=[top_sidebar, main, secondary]
            ),
            html.Div(
                id="bottom",
                className="bottom",
                children=[
                    # bottom_sidebar, bottom_figure
                ],
            ),
        ]
        return layout

    def _get_apply_callbacks(self):
        outputs = [
            Output(component_id="volcano-groupby", component_property="options"),
            Output(component_id="volcano-groupby", component_property="value"),
            Output(component_id="volcano-reference", component_property="options"),
            Output(component_id="volcano-reference", component_property="value"),
            Output(
                component_id="volcano-groupby-container", component_property="style"
            ),
            Output(
                component_id="volcano-reference-container", component_property="style"
            ),
        ]
        inputs = {
            "submit": Input(component_id="de-submit", component_property="n_clicks"),
        }
        states = {
            "groupby": State(component_id="de-groupby", component_property="value"),
        }
        for param in DE.de_params.values():
            states[param.key] = State(
                component_id=f"de-{param.key}", component_property="value"
            )

        return dict(output=outputs, inputs=inputs, state=states)

    def _get_plot_callbacks(self):
        outputs = [
            Output(component_id="de-volcano-plot", component_property="figure"),
            Output(component_id="de-secondary-plot", component_property="figure"),
        ]
        inputs = {
            "groupby": Input(
                component_id="volcano-groupby", component_property="value"
            ),
            "reference": Input(
                component_id="volcano-reference", component_property="value"
            ),
        }
        states = {}
        return dict(output=outputs, inputs=inputs, state=states)

    def _params_layout(self):
        cats = self.dataset.get_categoricals()
        divs = []

        groups = DE.get_groupby_options(dataset=self.dataset)
        refs = DE.get_reference_options(
            dataset=self.dataset, groupby=next(iter(groups), None)
        )
        # Volcano Group By Select
        divs.append(
            html.Div(
                id="volcano-groupby-container",
                children=[
                    html.Label("Volcano GroupBy"),
                    html.Div(
                        [
                            dcc.Dropdown(
                                options=groups,
                                value=next(iter(groups), None),
                                id="volcano-groupby",
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
                id="volcano-reference-container",
                children=[
                    html.Label("Volcano Reference"),
                    html.Div(
                        [
                            dcc.Dropdown(
                                options=refs,
                                value=next(iter(refs), None),
                                id="volcano-reference",
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

        layout = html.Div(children=divs)
        return layout

    def _get_on_click_callbacks(self):
        outputs = [Output("de-volcano-info", "children")]
        inputs = {
            "click_data": Input("de-volcano-plot", "clickData"),
        }
        return dict(output=outputs, inputs=inputs)
