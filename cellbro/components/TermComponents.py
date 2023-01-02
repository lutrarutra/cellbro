import dash
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, ctx, ALL

from ..util.DashAction import DashAction
from .DashComponent import DashComponent
from ..components.CID import CID

# class AddTermToList(DashAction):
#     def __init__(self, dataset, page_id):
#         super().__init__(dataset, page_id)

#     def setup_callbacks(self, app):

class HandleTermPopup(DashAction):
    def setup_callbacks(self, app):
        output = Output("termlist-popup", "is_open")

        inputs = dict(
            open=Input(dict(type="add-term-to-genelist", index=ALL), "n_clicks"),
            close=Input(f"termlist-popup-close", "n_clicks"),
            submit=Input("termlist-popup-new-genelist", "n_clicks"),
            name=State("termlist-genelist-name", "value"),
            genes=State("termlist-genes", "value")
        )

        @app.dash_app.callback(output=output, inputs=inputs)
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
                if name in self.dataset.adata.uns["genelists"].keys():
                    return True
                else:
                    self.dataset.adata.uns["genelists"][name] = genes
                
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
            genelists = self.dataset.get_genelists()
            #enabled = term_name in genelists
            return [title, term_name, genes, genes, genelists]


class AddGenesFromTermPopUp(DashComponent):
    def __init__(self, dataset):
        super().__init__("static", "static", "add_genes_from_term_popup")
        self.dataset = dataset
        self.actions.update(
            update_term_popup=UpdateTermPopup(CID(self.page_id, self.loc_class, "update_term_popup"), dataset),
            handle_popup=HandleTermPopup(CID(self.page_id, self.loc_class, "handle_popup"), dataset)
        )

    def create_layout(self):
        lists = self.dataset.get_genelists()

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


# def create_term_card(term_name, gene_set_name, genes):
#     if term_name is None:
#         return html.Div([
#             html.H4("Select Term by Clicking on a Point"),
#         ], style={
#             "display": "flex", "justify-content": "center",
#             "align-items": "center", "height": "100%", "width": "100%"
#         })

