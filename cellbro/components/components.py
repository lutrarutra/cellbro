from typing import Literal

import dash
from dash import html, dcc, Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import scout

from ..util.DashAction import DashAction
from .DashComponent import DashComponent
from .CID import LocClass

# continuous_colormaps = [
#     {"value" : "seismic", "label": "Seismic (for centered)"},
#     {"value" : "RdBu_r", "label": "B-W-R"},
#     {"value" : "viridis", "label": "Viridis"},
#     {"value" : "plasma", "label": "Plasma"},
#     {"value" : "inferno", "label": "Inferno"},
#     {"value" : "magma", "label": "Magma"},
#     {"value" : "cividis", "label": "Cividis"},
# ]


continuous_colormaps = list(scout.ply.get_continuous_colorscales().keys())
discrete_colormaps = list(scout.ply.get_discrete_colorscales().keys())


def create_colormap_selector(id, options, default=None):
    if default is None:
        default = list(options.keys())[0]

    return dcc.Dropdown(
        id=id,
        options=options,
        value=default,
        clearable=False,
    )

class FigureHeaderTab(DashComponent):
    def __init__(self, page_id: str, loc_class: LocClass, children, tab_label: str, id=None):
        super().__init__(page_id, loc_class, "fig_header_tab")
        self.children = children
        self.tab_label = tab_label
        self.id = id if id is not None else ""

    def create_layout(self):
        return dbc.Card(dbc.CardBody(self.children, className="param-row-stacked", id=self.id))

class FigureHeader(DashComponent):
    def __init__(self, page_id: str, loc_class: LocClass, tabs: list[FigureHeaderTab]):
        super().__init__(page_id, loc_class, "fig_header_tab")
        self.tabs = tabs

    def create_layout(self):
        if len(self.tabs) == 1:
            return self.tabs[0].create_layout()

        tabs = []
        for tab in self.tabs:
            tabs.append(dbc.Tab(tab.create_layout(), label=tab.tab_label))

        return dbc.Tabs(tabs, className="row-tabs")

class CollapseDiv(DashAction):
    def __init__(self, page_id, div_id, btn_id):
        super().__init__(dataset=None, page_id=page_id, loc_class="static")
        self.div_id = div_id
        self.btn_id = btn_id

    def setup_callbacks(self, app):
        @app.dash_app.callback(
            output=Output(self.div_id, "is_open"),
            inputs=[Input(self.btn_id, "n_clicks")],
            state=[State(self.div_id, "is_open")],
        )
        def _(n, is_open):
            if n:
                return not is_open
            return is_open

class CollapsibleDiv(DashComponent):
    def __init__(self, page_id, div_id, children, collapse_btn_id, collapsed=True):
        super().__init__(page_id)
        self.div_id = div_id
        self.children = children
        self.collapse_btn_id = collapse_btn_id
        self.children = children
        self.collapsed = collapsed
        self.actions = dict(
            collapse_div=CollapseDiv(self.page_id, self.div_id, self.collapse_btn_id)
        )

    def create_layout(self):
        return dbc.Collapse(
            children=self.children,
            id=self.div_id,
            is_open=not self.collapsed,
        )


def params_layout(params, id_prefix):
    divs = []
    for key, param in params.items():
        if param.type == list:
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
                                    options=param.allowed_values,
                                    value=param.default,
                                    id=f"{id_prefix}-{key}",
                                    clearable=False,
                                )
                            ],
                            className="param-select",
                        ),
                    ],
                    className="param-row-stacked",
                )
            )
        else:
            divs.append(
                html.Div(
                    [
                        html.Label(
                            param.name,
                            className="param-label",
                        ),
                        dcc.Input(
                            id=f"{id_prefix}-{key}",
                            type=param.input_type,
                            value=param.value,
                            step=param.step if param.step != None else 0.1,
                            className="param-input",
                        ),
                    ],
                    className="param-row",
                )
            )
            
    return divs


