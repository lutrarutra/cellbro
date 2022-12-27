import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

from ..util.DashFigure import DashFigure
from ..util.DashAction import PlotAction
from ..util import Components
from ..plots import DE

import scout

class PlotVolcano(PlotAction):
    def plot(self, groupby, reference):
        key = f"rank_genes_{groupby}"

        if reference is None:
            reference = list(self.dataset.adata.uns[key].keys())[0]

        fig = scout.ply.marker_volcano(
            self.dataset.adata.uns[key][reference], layout=DE.figure_layout
        )

        return fig

    def setup_callbacks(self, app):
        outputs = Output(component_id=f"{self.page_id_prefix}-{self.loc_class}-plot", component_property="figure")
        
        inputs = {
            "groupby": Input(
                component_id=f"{self.page_id_prefix}-{self.loc_class}-volcano-groupby", component_property="value"
            ),
            "reference": Input(
                component_id=f"{self.page_id_prefix}-{self.loc_class}-volcano-reference", component_property="value"
            ),
        }

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(groupby, reference):
            if groupby is None:
                raise PreventUpdate
            return self.plot(groupby, reference)

class DEVolcanoFig(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.actions.update(
            plot_projection=PlotVolcano(self.dataset, self.page_id_prefix, self.loc_class),
            click_action=Components.SelectGene(dataset=self.dataset, page_id_prefix=self.page_id_prefix, loc_class=self.loc_class)
        )

    def create_layout(self):
        select_gene_tab = Components.FigureParamTab(self.page_id_prefix, tab_label="Gene",
            id=f"{self.page_id_prefix}-{self.loc_class}-genecard", children=[
            Components.create_gene_card(self.page_id_prefix, self.loc_class, None, self.dataset)
        ])

        fig_header = Components.FigureParams(self.page_id_prefix, tabs=[select_gene_tab])

        figure = html.Div(
            children=[
                html.Div(fig_header.create_layout(), className="fig-header"),

                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id=f"{self.page_id_prefix}-{self.loc_class}-plot", className=f"{self.loc_class}-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    className=f"{self.loc_class}-body",
                ),
            ],
            className=f"{self.loc_class}",
        )
        return figure

    def get_sidebar_params(self) -> list:
        cats = self.dataset.get_categoric()
        divs = []

        groups = DE.get_groupby_options(dataset=self.dataset)
        refs = DE.get_reference_options(
            dataset=self.dataset, groupby=next(iter(groups), None)
        )
        # Volcano Group By Select
        divs.append(
            html.Div(
                children=[
                    html.Label("Volcano GroupBy"),
                    html.Div(
                        [
                            dcc.Dropdown(
                                options=groups,
                                value=next(iter(groups), None),
                                id=f"{self.page_id_prefix}-{self.loc_class}-volcano-groupby",
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
            html.Div([
                html.Label("Volcano Reference"),
                html.Div(
                    [
                        dcc.Dropdown(
                            options=refs,
                            value=next(iter(refs), None),
                            id=f"{self.page_id_prefix}-{self.loc_class}-volcano-reference",
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
