import plotly.graph_objects as go
import plotly.express as px
from dash import html, dcc, Input, Output, State
import dash_bootstrap_components as dbc


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