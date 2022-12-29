import dash
from dash import html, dcc, Input, Output, State, ctx, ALL
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from ..util.DashAction import DashAction
from .DashComponent import DashComponent

class SelectGene(DashAction):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)

    def apply(self, click_data):
        gene = click_data["points"][0]["hovertext"]
        element = create_gene_card(self.page_id_prefix, self.loc_class, gene, self.dataset)
        return element

    def setup_callbacks(self, app):
        outputs = Output(f"{self.page_id_prefix}-{self.loc_class}-genecard", "children")
        inputs = dict(
            click_data=Input(f"{self.page_id_prefix}-{self.loc_class}-plot", "clickData"),
        )

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(click_data):
            if click_data is None:
                raise PreventUpdate
            return self.apply(click_data)

class CreateGeneListAction(DashAction):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.selected_gene_cid = dict(id="selected-gene", page_id=ALL, loc_class=ALL)
        self.create_gene_list_btn_cid = dict(id="create-gene-list-button", page_id=ALL, loc_class=ALL)
        self.select_gene_list_cid = dict(id="gene-list-dropdown", page_id=ALL, loc_class=ALL)

    def create_new_list(self, name, selected_gene, available_list):
        if name is None or name == "":
            raise PreventUpdate
        if name in self.dataset.get_gene_lists():
            raise PreventUpdate

        assert name not in self.dataset.adata.uns["gene_lists"].keys()
        
        assert len(selected_gene) == 1
        selected_gene = next(iter(selected_gene), None)

        assert selected_gene is not None

        self.dataset.adata.uns["gene_lists"][name] = [selected_gene]

        return [False, dict(update=True), None, available_list]

    def setup_callbacks(self, app):
        output = [
            Output("genelist-popup", "is_open"),
            Output("genelist-store", "data"),
            Output("genelist-name", "value"),      # always returns None, just to clear it
            Output("genelist-existing-lists", "children"),
        ]
        inputs = dict(
            submit=Input(f"genelist-popup-submit", "n_clicks"),
            open=Input(self.create_gene_list_btn_cid, "n_clicks"),
            close=Input(f"genelist-popup-close", "n_clicks"),
        )
        state = dict(
            name=State(f"genelist-name", "value"),
            is_open=State(f"genelist-popup", "is_open"),
            selected_gene=State(self.selected_gene_cid, "children")
        )

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, open, close, name, is_open, selected_gene):
            # When a gene is selected new button is created, this prevents to go further without pressing it
            if all(v is None for v in open): raise PreventUpdate

            available_list = [dbc.ListGroupItem(name) for name in self.dataset.get_gene_lists()]

            if ctx.triggered_id == f"genelist-popup-close":
                return [False, dash.no_update, None, dash.no_update]

            if ctx.triggered_id == f"genelist-popup-submit":
                return self.create_new_list(name, selected_gene, available_list)

            return [True, dash.no_update, None, available_list]


class UpdateGeneList(DashAction):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.selected_gene_cid = dict(id="selected-gene", page_id=self.page_id_prefix, loc_class=self.loc_class)
        self.gene_list_dropdown_cid = dict(id="gene-list-dropdown", page_id=self.page_id_prefix, loc_class=self.loc_class)

    def setup_callbacks(self, app):
        output = [
            Output(dict(id="gene-list-dropdown", page_id=self.page_id_prefix, loc_class=self.loc_class), "options"),
            Output(dict(id="gene-list-dropdown", page_id=self.page_id_prefix, loc_class=self.loc_class), "value"),
        ]
        inputs = dict(
            genelist_store=Input(f"genelist-store", "data"),
            selected_list=Input(self.gene_list_dropdown_cid, "value"),
        )
        state = dict(
            selected_gene=State(self.selected_gene_cid, "children")
        )

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(genelist_store, selected_list, selected_gene):
            if ctx.triggered_id is None:
                raise PreventUpdate

            options = self.dataset.get_gene_lists()

            # if "update" in genelist_store.keys():
            #     selected_gene = genelist_store["selected_gene"]

            if ctx.triggered_id == self.gene_list_dropdown_cid:
                self.dataset.update_gene_lists(selected_gene, selected_list)

            value = self.dataset.get_gene_lists(selected_gene) if selected_gene is not None else []
            return options, value


class CreateGeneListPopup(DashComponent):
    def __init__(self, dataset):
        super().__init__("static")
        self.dataset = dataset
        self.actions.update(
            create_gene_list=CreateGeneListAction(self.dataset, self.page_id_prefix, loc_class="static"),
            # add_gene_to_list=AddGeneToList(self.dataset, self.page_id_prefix),
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
                    dbc.Button("Dummy", id=dict(id="create-gene-list-button", page_id=self.page_id_prefix, loc_class="dummy")),
                ], style={"display": "none"})
            ]),
            dbc.ModalFooter([
                dbc.Button("Create", id=f"genelist-popup-submit", n_clicks=0, color="primary"),
                dbc.Button("Close", id=f"genelist-popup-close", n_clicks=0, color="danger"),
            ])
        ])


def create_gene_card(page_id, loc_class, gene, dataset):
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
            html.H3(gene, id=dict(id="selected-gene", page_id=page_id, loc_class=loc_class)),
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
                    id=dict(id="gene-list-dropdown", page_id=page_id, loc_class=loc_class),
                    clearable=False,
                    placeholder="Select Gene List(s)",
                    multi=True,
                ),
            ], className="param-row-stacked", style={"width": "calc(100% - 100px)"}),
            html.Div([
                html.Label("New List"),
                html.Div([
                    dbc.Button("Create", id=dict(id="create-gene-list-button", page_id=page_id, loc_class=loc_class), color="primary")
                ]),
            ], className="param-row-stacked", style={"width": "100px"}),
        ], className="hover-body"),
    ], className="hover-container", style={"display": "flex"})

    return element
