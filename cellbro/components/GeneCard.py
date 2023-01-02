from collections import namedtuple

from dash import html, dcc, Output, Input, State, ctx
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from .DashComponent import DashComponent
from ..util.DashAction import DashAction
from .CID import CID

_genecards_query = "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
_google_scholar_query = "https://scholar.google.com/scholar?q="

class SelectGene(DashAction):
    RType = namedtuple(
        "RType",
        [
            "selected_gene",
            "genecards_link",
            "google_scholar_link",
            "genelist_options",
            "genelist_value",
            "placeholder_style",
            "card_style",
        ]
    )

    def __init__(self, cid: CID, dataset):
        super().__init__(cid, dataset)
        self.selected_gene_label_cid = dict(type="selected-gene", page_id=self.page_id, loc_class=self.loc_class.name)
        self.genelist_dropdown_cid = dict(type="select-genelist", page_id=self.page_id, loc_class=self.loc_class.name)

    def update_layout(self, gene):
        return self.RType(
            selected_gene=gene,
            genecards_link=f"{_genecards_query}{gene}",
            google_scholar_link=f"{_google_scholar_query}{gene}",
            genelist_options=self.dataset.get_genelists(),
            genelist_value=self.dataset.get_genelists(gene),
            placeholder_style=dict(display="none"),
            card_style=dict(display="flex", width="100%"),
        )

    def setup_callbacks(self, app):
        outputs = [
            Output(self.selected_gene_label_cid, "children"),
            Output(f"{self.page_id}-{self.loc_class}-genecards-link", "href"),
            Output(f"{self.page_id}-{self.loc_class}-google_scholar-link", "href"),
            Output(self.genelist_dropdown_cid, "options"),
            Output(self.genelist_dropdown_cid, "value"),
            Output(f"{self.page_id}-{self.loc_class}-placeholder", "style"),
            Output(f"{self.page_id}-{self.loc_class}-card", "style")
        ]
        inputs = dict(
            click_data=Input(dict(page_id=self.page_id, loc_class=self.loc_class.name, type="plot"), "clickData"),
            genelist_store=Input(f"genelist-store", "data"),
            selected_list=Input(self.genelist_dropdown_cid, "value"),
        )

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(click_data, genelist_store, selected_list):
            if click_data is None:
                raise PreventUpdate

            gene = click_data["points"][0]["hovertext"]
            if ctx.triggered_id == self.genelist_dropdown_cid:
                self.dataset.update_genelists(gene, selected_list)
            
            return self.update_layout(gene)


class GeneCard(DashComponent):
    def __init__(self, page_id, loc_class, dataset):
        super().__init__(page_id, loc_class, "gene_card")
        self.dataset = dataset
        self.actions.update(
            select_gene=SelectGene(self.cid, self.dataset)
        )

    def create_layout(self):
        gl_options = self.dataset.get_genelists()
        gl_elements = []
        for gl in gl_options:
            gl_elements.append({"label": gl, "value": gl})

        gl_chosen = []

        layout = html.Div([
            html.Div(
                id=f"{self.page_id}-{self.loc_class}-placeholder",
                children=html.H4("Select Gene by Clicking on a Point"),
                style={
                    "display": "flex", "justify-content": "center",
                    "align-items": "center", "height": "100%", "width": "100%"
                }
            ),
            html.Div(
                id=f"{self.page_id}-{self.loc_class}-card",
                style=dict(display="none"),
                children=[
                    html.Div([
                        html.Label("Gene:"),
                        html.H3("", id=dict(type="selected-gene", page_id=self.page_id, loc_class=self.loc_class.name)),
                    ], className="hover-header", style=dict(width="120px")),
                    html.Div([
                        html.Div([
                            html.A(
                                id=f"{self.page_id}-{self.loc_class}-genecards-link",
                                href=_genecards_query, role="button", target="_blank",
                                children=[html.Img(src="assets/logos/genecards_logo.png", style={"height": "20px"})]
                            ),
                            html.A(
                                id=f"{self.page_id}-{self.loc_class}-google_scholar-link",
                                href=_google_scholar_query, role="button", target="_blank",
                                children=[html.Img(src="assets/logos/google_scholar_logo.png", style={"height": "20px"})]
                            ),
                        ], className="hover-links"),
                        html.Div([
                            html.Label("Gene List(s)"),
                            dcc.Dropdown(
                                options=gl_elements,
                                value=gl_chosen,
                                id=dict(page_id=self.page_id, loc_class=self.loc_class.name, type="select-genelist"),
                                clearable=False,
                                placeholder="Select Gene List(s)",
                                multi=True,
                            ),
                        ], className="param-row-stacked", style={"width": "calc(100% - 100px)"}),
                        html.Div([
                            html.Label("New List"),
                            html.Div([
                                dbc.Button("Create", id=dict(page_id=self.page_id, loc_class=self.loc_class.name, type="create-genelist-btn"), color="primary")
                            ]),
                        ], className="param-row-stacked", style={"width": "100px"}),
                    ], className="hover-body"),
                ]
            )
        ], className="hover-container", style=dict(display="flex"))

        return layout
