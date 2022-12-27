from abc import ABC, abstractmethod
from typing import Literal

import dash
from dash import html, dcc, Input, Output, State, ctx, ALL
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import plotly
import plotly.express as px
import scanpy as sc

import scout

from cellbro.util.DashAction import DashAction

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

class SelectGene(DashAction):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix)
        self.loc_class = loc_class

    def apply(self, click_data):
        gene = click_data["points"][0]["hovertext"]
        element = create_gene_card(gene, self.dataset)
        return element

    def setup_callbacks(self, app):
        outputs = Output(f"{self.page_id_prefix}-{self.loc_class}-genecard", "children")
        inputs = {
            "click_data": Input(f"{self.page_id_prefix}-{self.loc_class}-plot", "clickData"),
        }

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(click_data):
            if click_data is None:
                raise PreventUpdate
            return self.apply(click_data)

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

class CreateGeneListAction(DashAction):
    def setup_callbacks(self, app):
        output = [
            Output(f"genelist-popup", "is_open"),
            Output(f"genelist-store", "data")
        ]
        inputs = dict(
            submit=Input(f"genelist-popup-submit", "n_clicks"),
            open=Input(dict(type="new-gene-list-button", index=ALL), "n_clicks"),
            close=Input(f"genelist-popup-close", "n_clicks"),
        )
        state = dict(
            name=State(f"genelist-name", "value"),
            is_open=State(f"genelist-popup", "is_open"),
            selected_gene=State(dict(type="selected-gene", index=ALL), "children")
        )
    
        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, open, close, name, is_open, selected_gene):
            if ctx.triggered_id is None:
                raise PreventUpdate

            if ctx.triggered_id == f"genelist-popup-close":
                return [False, dict()]

            if isinstance(ctx.triggered_id, dict):
                for v in open:
                    if v is not None:
                        return [True, dict()]

                return [False, dict()]

            if ctx.triggered_id == f"genelist-popup-submit":
                if name is None or name == "":
                    raise PreventUpdate
                if name in self.dataset.get_gene_lists():
                    raise PreventUpdate

                selected_gene = next(iter(selected_gene), None)
                self.dataset.adata.uns["gene_lists"][name] = [selected_gene] if selected_gene is not None else []

                return [False, dict(update=True, selected_gene=selected_gene)]

class NewGeneList(DashAction):
    def __init__(self, dataset, page_id_prefix):
        super().__init__(dataset, page_id_prefix)

    def setup_callbacks(self, app):
        output = [
            Output(dict(type="gene-list-dropdown", index=ALL), "options"),
            Output(dict(type="gene-list-dropdown", index=ALL), "value"),
        ]
        inputs = dict(
            genelist_store=Input(f"genelist-store", "data")
        )
        state = dict(
            dropdowns = State(dict(type="gene-list-dropdown", index=ALL), "id"),
        )
        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(genelist_store, dropdowns):
            if "update" not in genelist_store.keys():
                raise PreventUpdate

            options = [self.dataset.get_gene_lists()] * len(dropdowns)
            selected_gene = genelist_store["selected_gene"]
            value = self.dataset.get_gene_lists(selected_gene) if selected_gene is not None else []
            value = [value] * len(dropdowns)

            return options, value

class CreateGeneListPopup(DashComponent):
    def __init__(self, page_id_prefix, dataset):
        super().__init__(page_id_prefix)
        self.dataset = dataset
        self.actions.update(
            create_gene_list=CreateGeneListAction(self.dataset, self.page_id_prefix),
            new_gene_list=NewGeneList(self.dataset, self.page_id_prefix)
        )

    def create_layout(self):
        return dbc.Modal(id=f"genelist-popup", is_open=False, children=[
            dbc.ModalHeader("Create Gene List"),
            dbc.ModalBody([
                html.Div([
                    html.Label("Name:", className="param-label"),
                    html.Div([
                        dbc.Input(id=f"genelist-name", type="text", placeholder="Name"),
                    ], className="param-select")
                ], className="param-row-stacked"),
                # Dummy
                html.Div([
                    dbc.Button("Dummy", id=dict(type="new-gene-list-button", index=0)),
                ], style={"display": "none"})
            ]),
            dbc.ModalFooter([
                dbc.Button("Create", id=f"genelist-popup-submit", n_clicks=0, color="primary"),
                dbc.Button("Close", id=f"genelist-popup-close", n_clicks=0, color="danger"),
            ])
        ])


def create_gene_card(gene, dataset):
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


