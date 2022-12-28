import dash
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, ctx, ALL

from ..util.DashAction import DashAction
from . import components

# class AddTermToList(DashAction):
#     def __init__(self, dataset, page_id_prefix):
#         super().__init__(dataset, page_id_prefix)

#     def setup_callbacks(self, app):

class HandleTermPopup(DashAction):
    def setup_callbacks(self, app):
        output = Output("termlist-popup", "is_open")

        inputs = dict(
            open=Input(dict(type="add-term-to-gene-list", index=ALL), "n_clicks"),
            close=Input(f"termlist-popup-close", "n_clicks"),
            submit=Input("termlist-popup-new-genelist", "n_clicks"),
        )
        state = dict(
            name=State("termlist-genelist-name", "value"),
            genes=State("termlist-genes", "value")
        )

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(open, close, submit, name, genes):
            if ctx.triggered_id is None:
                return False

            if ctx.triggered_id == f"termlist-popup-close":
                return False

            if isinstance(ctx.triggered_id, dict):
                for v in open:
                    if v is not None:
                        return True
            
            if submit:
                if name in self.dataset.adata.uns["gene_lists"].keys():
                    return True
                else:
                    self.dataset.adata.uns["gene_lists"][name] = genes
                
            return False


class UpdateTermPopup(DashAction):
    def setup_callbacks(self, app):
        output = [
            Output("termlist-popup-term-title", "children"),
            Output("termlist-genelist-name", "value"),
            Output("termlist-genes", "options"),
            Output("termlist-genes", "value"),
            Output("termlist-existing-genelists", "options"),
            # Output("termlist-popup-new-genelist", "disabled")
        ]
        inputs = dict(
            store=Input("term-store", "data"),
        )
        @app.dash_app.callback(output=output, inputs=inputs) 
        def _(store):
            if store is None:
                raise PreventUpdate

            gene_set_name = store["gene_set_name"]
            genes = store["genes"]
            term_name = store["term_name"]

            title = f"{gene_set_name}: {term_name}"
            genelists = self.dataset.get_gene_lists()
            #enabled = term_name in genelists
            return [title, term_name, genes, genes, genelists]


class AddGenesFromTermPopUp(components.DashComponent):
    def __init__(self, page_id_prefix, dataset):
        super().__init__(page_id_prefix)
        self.dataset = dataset
        self.actions.update(
            update_term_popup=UpdateTermPopup(dataset, page_id_prefix),
            handle_popup=HandleTermPopup(dataset, page_id_prefix)
        )

    def create_layout(self):
        lists = self.dataset.get_gene_lists()

        return dbc.Modal(id=f"termlist-popup", is_open=False, children=[
            dbc.ModalHeader("Term: ", id="termlist-popup-term-title"),
            dbc.ModalBody([
                html.Div([
                    html.Label("Select Gene(s):", className="param-label"),
                    dcc.Dropdown(
                        options=[], value=None,
                        id="termlist-genes",
                        clearable=False,
                        placeholder="Select Gene(s)",
                        multi=True,
                    ),
                ], className="param-row-stacked"),

                html.Div([
                    html.Label("Create new List:", className="param-label"),
                    html.Div([
                        html.Div([
                            dbc.Input(
                                id="termlist-genelist-name", type="text", placeholder="Gene List Name",
                                value=""
                            ),
                            dbc.Button("Create", id="termlist-popup-new-genelist", n_clicks=None, color="primary")
                        ], style=dict(display="flex", gap="10px"))
                    ], className="param-select"),
                ], className="param-row-stacked"),

                html.Div([
                    html.Label("Existing Lists:", className="param-label"),
                    dcc.Dropdown(
                        options=lists, value=None,
                        id="termlist-existing-genelists",
                        clearable=False, placeholder="Select Gene List(s)",
                        multi=True
                    )
                ], className="param-row-stacked"),
            ]),
            dbc.ModalFooter([
                # dbc.Button("Create", id=f"termlist-popup-submit", n_clicks=0, color="primary"),
                dbc.Button("Close", id=f"termlist-popup-close", n_clicks=None, color="danger"),
            ])
        ])

class SelectTerm(DashAction):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix)
        self.loc_class = loc_class

    def apply(self, click_data, groupby, reference):
        term = click_data["points"][0]["hovertext"]
        df = self.dataset.adata.uns["gsea"][groupby][reference]
        gene_set_name = df.index.name
        genes = df[df["Term"] == term]["lead_genes"].values.tolist()[0]
        element = create_term_card(term, gene_set_name, genes)

        return element, dict(term_name=term, gene_set_name=gene_set_name, genes=genes)

    def setup_callbacks(self, app):
        outputs = [
            Output(f"{self.page_id_prefix}-{self.loc_class}-termcard", "children"),
            Output("term-store", "data")
        ]
        inputs = dict(
            click_data=Input(f"{self.page_id_prefix}-{self.loc_class}-plot", "clickData"),
        )
        state = dict(
            groupby=State(component_id=f"{self.page_id_prefix}-{self.loc_class}-groupby", component_property="value"),
            reference=State(component_id=f"{self.page_id_prefix}-{self.loc_class}-reference", component_property="value"),
        )

        @app.dash_app.callback(output=outputs, inputs=inputs, state=state)
        def _(click_data, groupby, reference):
            if click_data is None:
                raise PreventUpdate
            return self.apply(click_data, groupby, reference)


def create_term_card(term_name, gene_set_name, genes):
    if term_name is None:
        return html.Div([
            html.H4("Select Term by Clicking on a Point"),
        ], style={
            "display": "flex", "justify-content": "center",
            "align-items": "center", "height": "100%", "width": "100%"
        })

    element = html.Div([
        html.Div([
            html.Label(f"Term ({gene_set_name}):"),
            html.H3(term_name, id=dict(type="selected-term", index=0)),
        ], className="hover-header", style=dict(width="200px")),
        html.Div([
            html.Div([
                # html.A(
                #     href=f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}", role="button", target="_blank",
                #     children=[html.Img(src="assets/logos/genecards_logo.png", style={"height": "20px"})]
                # ),
                html.A(
                    href=f"https://scholar.google.com/scholar?q={term_name}", role="button", target="_blank",
                    children=[html.Img(src="assets/logos/google_scholar_logo.png", style={"height": "20px"})]
                ),
            ], className="hover-links"),
            html.Div([
                html.Label("Genes"),
                dcc.Dropdown(
                    options=genes,
                    value=genes,
                    id=dict(type="term-gene-list-dropdown", index=0),
                    clearable=False,
                    placeholder="Select Gene List(s)",
                    multi=True,
                ),
            ], className="param-row-stacked long-dropdown", style={"width": "calc(100% - 120px)"}),
            html.Div([
                html.Label("Add to Gene List"),
                html.Div([
                    dbc.Button("Add", id=dict(type="add-term-to-gene-list", index=0), color="primary", style=dict(width="100%"))
                ]),
            ], className="param-row-stacked", style={"width": "120px"}),
        ], className="hover-body"),
    ], className="hover-container", style={"display": "flex"})

    return element
