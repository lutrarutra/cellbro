from collections import namedtuple

from dash import html, dcc, Output, Input, State, ctx
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from .DashComponent import DashComponent
from ..util.DashAction import DashAction
from .CID import CID

_scholar_url_prefix = "https://scholar.google.com/scholar?q="


class SelectTerm(DashAction):
    RType = namedtuple(
        "RType",
        [
            "selected_term",
            "term_link",
            "term_genes",
            "geneset_name",
            "placeholder_style",
            "card_style",
        ]
    )

    def __init__(
        self, cid: CID, dataset,
        select_groupby_cid: CID,
        select_reference_cid: CID,
        select_library_cid: CID,
        selected_geneset_cid: CID = None,
        selected_term_cid: CID = None,
        term_genes_cid: CID = None,
        term_link_cid: CID = None,
    ):
        super().__init__(cid, dataset)
        self.plot_cid = CID(self.parent_cid.page_id, self.parent_cid.loc_class, "plot")
        self.select_groupby_cid = select_groupby_cid
        self.select_reference_cid = select_reference_cid
        self.select_library_cid = select_library_cid

        if selected_geneset_cid is not None:
            self.selected_geneset_cid = selected_geneset_cid
        else:
            self.selected_geneset_cid = CID(self.parent_cid.page_id, self.parent_cid.loc_class, "selected-geneset")

        if selected_term_cid is not None:
            self.selected_term_cid = selected_term_cid
        else:
            self.selected_term_cid = CID(self.parent_cid.page_id, self.parent_cid.loc_class, "selected-term")

        if term_genes_cid is not None:
            self.term_genes_cid = term_genes_cid
        else:
            self.term_genes_cid = CID(self.parent_cid.page_id, self.parent_cid.loc_class, "term-genes")

        if term_link_cid is not None:
            self.term_link_cid = term_link_cid
        else:
            self.term_link_cid = CID(self.parent_cid.page_id, self.parent_cid.loc_class, "term-link")

    def apply(self, click_data, groupby, reference, library):
        term = click_data["points"][0]["hovertext"]
        print(groupby)
        print(reference)
        print(library)
        print(self.dataset.adata.uns["gsea"].keys())
        print(self.dataset.adata.uns["gsea"][groupby].keys())
        df = self.dataset.adata.uns["gsea"][groupby][reference][library]
        geneset_name = df.index.name
        term_genes = df[df["Term"] == term]["lead_genes"].values.tolist()[0]

        return self.RType(
            selected_term=term,
            term_link=f"{_scholar_url_prefix}{term}",
            term_genes=term_genes,
            geneset_name=geneset_name,
            placeholder_style=dict(display="none"),
            card_style=dict(display="flex", width="100%"),
        )

    def setup_callbacks(self, app):
        outputs = [
            Output(self.selected_term_cid.to_dict(), "children"),
            Output(self.term_link_cid.to_dict(), "href"),
            Output(self.term_genes_cid.to_dict(), "options"),
            Output(self.selected_geneset_cid.to_dict(), "children"),
            Output(f"{self.page_id}-{self.loc_class}-placeholder", "style"),
            Output(f"{self.page_id}-{self.loc_class}-card", "style"),
        ]
        inputs = dict(
            click_data=Input(self.plot_cid.to_dict(), "clickData"),
            groupby=State(self.select_groupby_cid.to_dict(), "value"),
            reference=State(self.select_reference_cid.to_dict(), "value"),
            library=State(self.select_library_cid.to_dict(), "value"),
        )

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(click_data, groupby, reference, library):
            if click_data is None:
                raise PreventUpdate

            return self.apply(click_data, groupby, reference, library)

class TermCard(DashComponent):
    def __init__(self, page_id, loc_class, dataset, select_groupby_cid, select_reference_cid, select_library_cid):
        super().__init__(page_id, loc_class, "term_card")
        self.dataset = dataset
        self.actions.update(
            select_term=SelectTerm(
                self.cid, self.dataset,
                select_groupby_cid, select_reference_cid, 
                select_library_cid
            )
        )

    def create_layout(self):
        element = html.Div([
            html.Div(
                id=f"{self.page_id}-{self.loc_class}-placeholder",
                children=html.H4("Select Term by Clicking on a Point"),
                style={
                    "display": "flex", "justify-content": "center",
                    "align-items": "center", "height": "100%", "width": "100%"
                }
            ),
            html.Div(
                id=f"{self.page_id}-{self.loc_class}-card", style=dict(display="none"), children=[
                    html.Div([
                        html.Label(f"Term (DUMMY)", id=dict(type="selected-geneset", page_id=self.page_id, loc_class=self.loc_class.name)),
                        html.H5("", id=dict(type="selected-term", page_id=self.page_id, loc_class=self.loc_class.name), style={"font-size": "15px"}),
                    ], className="hover-header", style=dict(width="200px")),
                    html.Div([
                        html.Div([
                            html.A(
                                href=_scholar_url_prefix, id=dict(type="term-link", page_id=self.page_id, loc_class=self.loc_class.name), role="button", target="_blank",
                                children=[html.Img(src="assets/logos/google_scholar_logo.png", style={"height": "20px"})]
                            ),
                        ], className="hover-links"),
                        html.Div([
                            html.Label("Genes"),
                            dcc.Dropdown(
                                options=[],
                                id=dict(type="term-genes", page_id=self.page_id, loc_class=self.loc_class.name),
                                clearable=False, placeholder="Select Gene List(s)",
                                multi=True,
                            ),
                        ], className="param-row-stacked long-dropdown", style={"width": "calc(100% - 120px)"}),
                        html.Div([
                            html.Label("Add to Gene List"),
                            html.Div([
                                dbc.Button("Add", id=dict(type="add-term-to-genelist", index=0), color="primary", style=dict(width="100%"))
                            ]),
                        ], className="param-row-stacked", style={"width": "120px"}),
                    ], className="hover-body"),
                ]
            )
        ], className="hover-container", style=dict(display="flex"))
        

        return element
