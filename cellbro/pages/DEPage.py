import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

import scout

from ..util.DashAction import DashAction
from ..components.DashPage import DashPage
from ..components.CID import CID
from ..components.Sidebar import Sidebar
from ..components import components
from ..plots.DE.DEHistogram import DEHistogram

from ..plots import DE

class ApplyDE(DashAction):
    def apply(self, groupby, reference, params):
        scout.tl.rank_marker_genes(self.dataset.adata, groupby=groupby, reference=reference, **params)

        return dict(update=True, groupby=groupby)

    def setup_callbacks(self, app):
        output = Output("de-store", "data")

        inputs = dict(
            submit=Input(f"{self.page_id}-main-sidebar-apply_btn", "n_clicks"),
            groupby=State(f"{self.page_id}-groupby", "value"),
            reference=State(f"{self.page_id}-reference", "value")
        )

        for param in DE.de_tools.de_params.values():
            inputs[param.key] = State(f"de-{param.key}", "value")

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(submit, groupby, reference, **kwargs):
            if submit is None:
                raise PreventUpdate

            return self.apply(groupby, reference, params=kwargs)

        # Update References
        output = Output(f"{self.page_id}-reference", "options")
        inputs = dict(
            groupby=Input(f"{self.page_id}-groupby", "value"),
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
                        dcc.Dropdown(options=cats, value=cats[0], id=f"{self.page_id}-groupby", clearable=False)
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
                            options=["rest"], value="rest", id=f"{self.page_id}-reference", clearable=False
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
                                id=f"{self.page_id}-{key}", value=param.default, options=param.allowed_values,
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
        self.components.update(
            main_figure=DE.DEVolcano(dataset=self.dataset, page_id=self.page_id, loc_class="main"),
        )
        self.components.update(
            secondary_figure=DEHistogram(
                self.dataset, page_id=self.page_id, loc_class="secondary",
                select_groupby_cid=self.components["main_figure"].children["select_groupby"].cid,
                select_reference_cid=self.components["main_figure"].children["select_reference"].cid,
            )
        )
        self.actions.update(
            de_apply=ApplyDE(CID(self.page_id, "static", "de_apply"), self.dataset),
        )

    def create_layout(self):
        self.components["top_sidebar"] = Sidebar(
            page_id=self.page_id, loc_class="main",
            title="Differential Expression Settings", create_btn=True,
            params_children=self.actions["de_apply"].get_sidebar_params(),
        )

        self.components["bot_sidebar"] = Sidebar(
            page_id=self.page_id, loc_class="bottom",
            title="Empty", create_btn=False,
            params_children=[],
        )

        layout = [
            html.Div(
                className="top", children=[
                    self.components["top_sidebar"].create_layout(), 
                    self.components["main_figure"].create_layout(),
                    self.components["secondary_figure"].create_layout()
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
