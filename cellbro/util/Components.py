from abc import ABC, abstractmethod
from typing import Literal

import dash
from dash import html, dcc, Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import scout

from .DashAction import DashAction

continuous_colormaps = [
    {"value" : "seismic", "label": "Seismic (for centered)"},
    {"value" : "RdBu_r", "label": "B-W-R"},
    {"value" : "viridis", "label": "Viridis"},
    {"value" : "plasma", "label": "Plasma"},
    {"value" : "inferno", "label": "Inferno"},
    {"value" : "magma", "label": "Magma"},
    {"value" : "cividis", "label": "Cividis"},
]

discrete_colormaps = [
    {"value" : value, "label" : value.replace("_", " ").title()} for value in scout.ply._discrete_cmap_mapping.keys()
]


class DashComponent(ABC):
    def __init__(self, page_id_prefix):
        self.actions = {}
        self.page_id_prefix = page_id_prefix

    @abstractmethod
    def create_layout(self):
        ...

    def setup_callbacks(self, app):
        for action in self.actions.values():
            action.setup_callbacks(app)


def create_colormap_selector(id, options, default=None):
    if default is None:
        default = list(options.keys())[0]

    return dcc.Dropdown(
        id=id,
        options=options,
        value=default,
        clearable=False,
    )

class FigureParamTab(DashComponent):
    def __init__(self, page_id_prefix, children, tab_label, id=None):
        super().__init__(page_id_prefix)
        self.children = children
        self.tab_label = tab_label
        self.id = id if id is not None else ""

    def create_layout(self):
        return dbc.Card(dbc.CardBody(self.children, className="param-row-stacked", id=self.id))

class FigureParams(DashComponent):
    def __init__(self, page_id_prefix, tabs: list[FigureParamTab]):
        super().__init__(page_id_prefix)
        self.tabs = tabs

    def create_layout(self):
        if len(self.tabs) == 1:
            return self.tabs[0].create_layout()

        tabs = []
        for tab in self.tabs:
            tabs.append(dbc.Tab(tab.create_layout(), label=tab.tab_label))

        return dbc.Tabs(tabs, className="row-tabs")

class HideSidebar(DashAction):
    def __init__(self, page_id_prefix, id, btn_id, side: Literal["left", "right"]):
        super().__init__(dataset=None, page_id_prefix=page_id_prefix)
        self.id = id
        self.btn_id = btn_id
        self.side = side

    def setup_callbacks(self, app):
        @app.dash_app.callback(
            output=Output(self.id, "className"),
            inputs=[Input(self.btn_id, "value")],
            state=[State(self.id, "className")],
        )
        def _(value, cls):
            if value != None:
                if value:
                    return cls + f" {self.side}-open"
                else:
                    return cls.replace(f"{self.side}-open", "")

            raise PreventUpdate

class Sidebar(DashComponent):
    def __init__(
        self, page_id_prefix, title, params_children, apply_btn_id,
        side: Literal["left", "right"], row: Literal["top", "bot"], btn_text = "Apply"
    ):
        super().__init__(page_id_prefix)
        self.id = f"{page_id_prefix}-{side}-{row}-sidebar"
        self.apply_btn_id = apply_btn_id
        self.title = title
        self.params_children = params_children
        self.row = row
        self.side = side
        self.btn_text = btn_text
        self.actions = dict(
            hide_sidebar=HideSidebar(
                page_id_prefix=self.page_id_prefix, id=self.id,
                btn_id=f"{self.side}-sidebar-btn", side=self.side
            ),
        )

    def create_layout(self):
        if self.apply_btn_id is not None:
            btn_container = [
                dbc.Button(
                    self.btn_text,
                    color="primary",
                    className="mr-1",
                    id=self.apply_btn_id,
                ),
            ]
        else: btn_container = []
        
        return html.Div(
            children=[
                html.Div(
                    [
                        html.H3(self.title),
                    ],
                    className="sidebar-header",
                ),
                dcc.Loading(
                    type="circle", parent_className=f"sidebar-container",
                    children=[
                        html.Div(
                            children=self.params_children,
                            className=f"sidebar-parameters",
                        ),
                        html.Div(
                            children=btn_container,
                            className=f"sidebar-footer",
                        ),
                    ],
                ),
            ],
            className=f"{self.side}-sidebar sidebar {self.row}-sidebar", id=self.id
        )

class CollapseDiv(DashAction):
    def __init__(self, page_id_prefix, div_id, btn_id):
        super().__init__(dataset=None, page_id_prefix=page_id_prefix)
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
    def __init__(self, page_id_prefix, div_id, children, collapse_btn_id, collapsed=True):
        super().__init__(page_id_prefix)
        self.div_id = div_id
        self.children = children
        self.collapse_btn_id = collapse_btn_id
        self.children = children
        self.collapsed = collapsed
        self.actions = dict(
            collapse_div=CollapseDiv(self.page_id_prefix, self.div_id, self.collapse_btn_id)
        )

    def create_layout(self):
        return dbc.Collapse(
            children=self.children,
            id=self.div_id,
            is_open=not self.collapsed,
        )


def create_gene_card(gene, dataset):
    if gene is None:
        return html.Div([
            html.H4("Select Gene by Clicking on a Point"),
        ], style={
            "display": "flex", "justify-content": "center",
            "align-items": "center", "height": "100%", "width": "100%"
        })

    gl_options = dataset.get_gene_lists()
    gl_elements = []
    for gl in gl_options:
        gl_elements.append({"label": gl, "value": gl})
    
    gl_chosen = dataset.get_gene_lists(gene=gene)

    element = html.Div([
        html.Div([
            html.Label("Gene:"),
            html.H3(gene, id=dict(type="selected-gene", index=0)),
        ], className="hover-header"),
        html.Div([
            html.Div([
                html.A(
                    href=f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}", role="button", target="_blank",
                    children=[html.Img(src="assets/logos/genecards_logo.png", style={"height": "20px"})]
                ),
                html.A(
                    href=f"https://scholar.google.com/scholar?q={gene}", role="button", target="_blank",
                    children=[html.Img(src="assets/logos/google_scholar_logo.png", style={"height": "20px"})]
                ),
            ], className="hover-links"),
            html.Div([
                html.Label("Gene List(s)"),
                dcc.Dropdown(
                    options=gl_elements,
                    value=gl_chosen,
                    id=dict(type="gene-list-dropdown", index=0),
                    clearable=False,
                    placeholder="Select Gene List(s)",
                    multi=True,
                ),
            ], className="param-row-stacked", style={"width": "calc(100% - 100px)"}),
            html.Div([
                html.Label("New List"),
                html.Div([
                    dbc.Button("Create", id=dict(type="new-gene-list-button", index=0), color="primary")
                ]),
            ], className="param-row-stacked", style={"width": "100px"}),
        ], className="hover-body"),
    ], className="hover-container",
    style={"display": "none" if gene == None else "flex"})

    return element


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


