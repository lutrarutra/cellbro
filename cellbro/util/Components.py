import plotly.graph_objects as go
import plotly.express as px
from dash import html, dcc, Input, Output, State
import dash_bootstrap_components as dbc

from cellbro.util.DashAction import DashAction


class CollapseDiv(DashAction):
    def __init__(self, id, btn_id, children, collapsed=True):
        super().__init__(dataset=None)
        self.id = id
        self.btn_id = btn_id
        self.children = children
        self.collapsed = collapsed
        self.layout = self.create_layout()

    def create_layout(self):
        return dbc.Collapse(
            children=self.children,
            id=self.id,
            is_open=not self.collapsed,
        )

    def apply(self, params):
        pass

    def setup_callbacks(self, dash_app):
        @dash_app.callback(
            output=Output(self.id, "is_open"),
            inputs=[Input(self.btn_id, "n_clicks")],
            state=[State(self.id, "is_open")],
        )
        def _(n, is_open):
            if n:
                return not is_open
            return is_open


def create_gene_card(gene, dataset):
    gl_options = dataset.get_gene_lists()
    gl_elements = []
    for gl in gl_options:
        gl_elements.append({"label": gl, "value": gl})
    
    gl_chosen = dataset.get_gene_lists(gene=gene)

    element = html.Div([
        html.Div([
            html.Label("Gene:", className="hover-title"),
            html.H3(gene, id="selected-gene", className="hover-title"),
            html.Div([
                html.A(
                    href=f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}", role="button", target="_blank",
                    children=[html.Img(src="assets/logos/genecards_logo.png", style={"height": "20px"})]
                ),
                html.A(
                    href=f"https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q={gene}", role="button", target="_blank",
                    children=[html.Img(src="assets/logos/google_scholar_logo.png", style={"height": "20px"})]
                ),
            ], className="hover-row"),
        ], className="hover-header"),
        html.Div([
            html.Div([
                html.Label("Gene List(s)"),
                dcc.Dropdown(
                    options=gl_elements,
                    value=gl_chosen,
                    id="gene-list-dropdown",
                    clearable=False,
                    placeholder="Select Gene List(s)",
                    multi=True,
                    style={"flex": "1"},
                ),
            ], className="hover-col"),
            html.Div([
                html.Label("Create New Gene List"),
                html.Div([
                    dbc.Input(placeholder="Gene List Name", id="new-gene-list-input", type="text"),
                    dbc.Button("Create", id="new-gene-list-button", color="primary", className="mr-1")
                ], className="hover-row"),
            ], className="hover-col"),
        ], className="hover-info")
    ], className="hover-container", style={"display": "none" if gene == None else "flex"})
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


