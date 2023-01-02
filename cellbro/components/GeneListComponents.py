from collections import namedtuple

import dash
from dash import html, dcc, Input, Output, State, ctx, ALL
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from ..util.DashAction import DashAction
from .DashComponent import DashComponent
from ..components.CID import CID

class CreateGeneListAction(DashAction):
    RType = namedtuple("RType", ["is_open", "genelist_store", "genelist_name", "existing_lists"])

    def __init__(self, dataset, page_id, loc_class):
        super().__init__(CID(page_id, loc_class, "create_genelist"), dataset)
        self.selected_gene_cid = dict(type="selected-gene", page_id=ALL, loc_class=ALL)
        self.create_genelist_btn_cid = dict(type="create-genelist-btn", page_id=ALL, loc_class=ALL)
        self.select_genelist_cid = dict(type="select-genelist", page_id=ALL, loc_class=ALL)

    def setup_callbacks(self, app):
        output = [
            Output("genelist-popup", "is_open"),
            Output("genelist-store", "data"),
            Output("genelist-name", "value"),      # always returns None, just to clear it
            Output("genelist-existing-lists", "children"),
        ]
        inputs = dict(
            submit=Input(f"genelist-popup-submit", "n_clicks"),
            open=Input(self.create_genelist_btn_cid, "n_clicks"),
            close=Input(f"genelist-popup-close", "n_clicks"),
            name=State(f"genelist-name", "value"),
            is_open=State(f"genelist-popup", "is_open"),
            selected_gene=State(self.selected_gene_cid, "children")
        )

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(submit, open, close, name, is_open, selected_gene):
            # When a gene is selected new button is created, this prevents to go further without pressing it
            if all(v is None for v in open): raise PreventUpdate

            available_list = [dbc.ListGroupItem(name) for name in self.dataset.get_genelists()]

            if ctx.triggered_id == f"genelist-popup-close":
                return self.RType(
                    is_open=False,
                    genelist_store=dash.no_update,
                    genelist_name=None,
                    existing_lists=dash.no_update
                )

            if ctx.triggered_id == f"genelist-popup-submit":
                return self.create_new_list(name, selected_gene, available_list)

            return self.RType(
                is_open=True,
                genelist_store=dash.no_update,
                genelist_name=None,
                existing_lists=available_list
            )

    def create_new_list(self, name, selected_gene, available_list):
        if name is None or name == "":
            raise PreventUpdate
        if name in self.dataset.get_genelists():
            raise PreventUpdate

        assert name not in self.dataset.adata.uns["genelists"].keys()

        assert len(selected_gene) == 1
        selected_gene = next(iter(selected_gene), None)

        assert selected_gene is not None

        self.dataset.adata.uns["genelists"][name] = [selected_gene]

        return self.RType(
            is_open=False,
            genelist_store=dict(update=True),
            genelist_name=None,
            existing_lists=available_list
        )


class CreateGeneListPopup(DashComponent):
    def __init__(self, dataset):
        super().__init__("static", "static", "genelist_popup")
        self.dataset = dataset
        self.actions.update(
            create_genelist=CreateGeneListAction(self.dataset, self.page_id, loc_class="static"),
            # add_gene_to_list=AddGeneToList(self.dataset, self.page_id),
        )

    def create_layout(self):
        lists = self.dataset.get_genelists()
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
                    dbc.Button("Dummy", id=dict(page_id="static", loc_class="dummy", type="create-genelist-btn")),
                ], style={"display": "none"})
            ]),
            dbc.ModalFooter([
                dbc.Button("Create", id=f"genelist-popup-submit", n_clicks=0, color="primary"),
                dbc.Button("Close", id=f"genelist-popup-close", n_clicks=0, color="danger"),
            ])
        ])
