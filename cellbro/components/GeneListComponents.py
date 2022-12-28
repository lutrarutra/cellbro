import dash
from dash import html, dcc, Input, Output, State, ctx, ALL
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from ..util.DashAction import DashAction
from . import components

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

class CreateGeneListAction(DashAction):
    def setup_callbacks(self, app):
        output = [
            Output(f"genelist-popup", "is_open"),
            Output(f"genelist-store", "data"),
            Output(f"genelist-name", "value"),
            Output("genelist-existing-lists", "children")
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
            available_list = [dbc.ListGroupItem(name) for name in self.dataset.get_gene_lists()]

            if ctx.triggered_id is None:
                raise PreventUpdate

            if ctx.triggered_id == f"genelist-popup-close":
                return [False, dict(), None, available_list]

            if isinstance(ctx.triggered_id, dict):
                for v in open:
                    if v is not None:
                        return [True, dict(), None, available_list]

                return [False, dict(), None, available_list]

            if ctx.triggered_id == f"genelist-popup-submit":
                if name is None or name == "":
                    raise PreventUpdate
                if name in self.dataset.get_gene_lists():
                    raise PreventUpdate

                selected_gene = next(iter(selected_gene), None)
                self.dataset.adata.uns["gene_lists"][name] = [selected_gene] if selected_gene is not None else []

                return [False, dict(update=True, selected_gene=selected_gene), None, available_list]


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
            dropdowns=State(dict(type="gene-list-dropdown", index=ALL), "id"),
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


class CreateGeneListPopup(components.DashComponent):
    def __init__(self, page_id_prefix, dataset):
        super().__init__(page_id_prefix)
        self.dataset = dataset
        self.actions.update(
            create_gene_list=CreateGeneListAction(self.dataset, self.page_id_prefix),
            new_gene_list=NewGeneList(self.dataset, self.page_id_prefix)
        )

    def create_layout(self):
        lists = self.dataset.get_gene_lists()
        return dbc.Modal(id=f"genelist-popup", is_open=False, children=[
            dbc.ModalHeader("Create Gene List"),
            dbc.ModalBody([
                html.Div([
                    html.Label("Name:", className="param-label"),
                    html.Div([
                        dbc.Input(id=f"genelist-name", type="text", placeholder="Gene List Name"),
                    ], className="param-select"),
                ], className="param-row-stacked"),

                html.Div([
                    html.Label("Existing Lists:", className="param-label"),
                    dbc.ListGroup(id="genelist-existing-lists",
                        children=[dbc.ListGroupItem(name) for name in lists],
                        style={"max-height": "200px", "overflow-y": "scroll"}
                    ),
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
        ], className="hover-header", style=dict(width="120px")),
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
    ], className="hover-container", style={"display": "flex"})

    return element
