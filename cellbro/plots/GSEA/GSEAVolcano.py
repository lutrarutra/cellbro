import dash
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, ctx

import gseapy

from ...components.DashPlot import DashPlot
from ...util.DashAction import DashAction
from ...components import components
from ...components.CID import CID
from ...components.DropDown import DropDown
from ...plots.DE import de_tools
from ...components.TermCard import TermCard

import scout

default_layout = dict(
    paper_bgcolor="white",
    plot_bgcolor="white",
    margin=dict(t=5, b=5, l=5, r=5),
)

# class ApplyGSEA(DashAction):
#     def __init__(
#         self, parent_cid, dataset,
#         select_gsea_groupby_cid, select_gsea_reference_cid,
#         select_geneset_cid,
#     ):
#         super().__init__(parent_cid, dataset)
#         self.select_gsea_groupby_cid = select_gsea_groupby_cid
#         self.select_gsea_reference_cid = select_gsea_reference_cid
#         self.select_geneset_cid = select_geneset_cid

#     def apply(self, groupby, reference, gene_set):
#         # gene_score_df = self.dataset.adata.uns[f"rank_genes_{groupby}"][reference]
#         gene_score_df = self.dataset.adata.uns["de"][groupby][reference]
#         res = scout.tl.GSEA(gene_score_df, score_of_interest="gene_score", gene_set=gene_set)
#         self.dataset.add_gsea_result(res, groupby, reference)

#     def setup_callbacks(self, app):
#         output = Output("gsea-store", "data")
#         inputs = dict(
#             submit=Input(f"{self.page_id}-{self.loc_class}-sidebar-apply_btn", "n_clicks"),
#             gsea_groupby=State(self.select_gsea_groupby_cid.to_dict(), "value"),
#             gsea_reference=State(self.select_gsea_reference_cid.to_dict(), "value"),
#             gsea_geneset=State(self.select_geneset_cid.to_dict(), "value"),
#         )

#         @app.dash_app.callback(output=output, inputs=inputs, prevent_initial_call=True)
#         def _(submit, gsea_groupby, gsea_reference, gsea_geneset):
#             if submit is None:
#                 raise PreventUpdate
                
#             if gsea_groupby is None or gsea_reference is None or gsea_geneset is None:
#                 raise PreventUpdate

#             self.apply(gsea_groupby, gsea_reference, gsea_geneset)
#             return dict(update=True)


class Plot(DashAction):
    def __init__(
        self, parent_cid, dataset,
        select_groupby_cid, select_reference_cid,
        select_geneset_cid,
    ):
        super().__init__(parent_cid, dataset)
        self.select_groupby_cid = select_groupby_cid
        self.select_reference_cid = select_reference_cid
        self.select_geneset_cid = select_geneset_cid

    def plot(self, groupby, reference, geneset):
        res = self.dataset.adata.uns[f"gsea"][groupby][reference][geneset]
        return scout.ply.gsea_volcano(res, layout=default_layout)

    def setup_callbacks(self, app):
        output = Output(self.parent_cid.to_dict(), "figure")
        inputs = dict(
            groupby=Input(self.select_groupby_cid.to_dict(), "value"),
            reference=Input(self.select_reference_cid.to_dict(), "value"),
            geneset=Input(self.select_geneset_cid.to_dict(), "value"),
        )
        @app.dash_app.callback(output=output, inputs=inputs) 
        def _(groupby, reference, geneset):
            if groupby is None or reference is None:
                raise PreventUpdate

            return self.plot(groupby, reference, geneset)

class GSEAVolcano(DashPlot):
    def __init__(self, dataset, page_id, loc_class):
        super().__init__(dataset, page_id, loc_class)

        def get_groupby_options():
            return sorted(list(self.dataset.adata.uns["gsea"].keys()))

        groupby_options = get_groupby_options()

        self.children["select_groupby"] = DropDown(
            cid=CID(page_id, loc_class, "select-groupby"),
            options=groupby_options, default=next(iter(groupby_options), None),
            options_callback=lambda: get_groupby_options(),
            update_store_id="gsea-store"
        )

        def get_reference_options(groupby):
            ref_options = sorted(list(self.dataset.adata.uns["gsea"][groupby].keys()))
            return ref_options

        if len(groupby_options) > 0:
            ref_options = get_reference_options(groupby_options[0])
        else:
            ref_options = []

        self.children["select_reference"] = DropDown(
            cid=CID(page_id, loc_class, "select-reference"),
            options=ref_options, default=next(iter(ref_options), None),
            options_callback=lambda x: get_reference_options(x),
            master_cid=self.children["select_groupby"].cid,
        )

        def get_geneset_options(groupby, target):
            ref_options = sorted(list(self.dataset.adata.uns["gsea"][groupby][target].keys()))
            return ref_options

        if len(ref_options) > 0:
            geneset_options = get_geneset_options(groupby_options[0], ref_options[0])
        else:
            geneset_options = []

        self.children["select_library"] = DropDown(
            cid=CID(page_id, loc_class, "select-library"),
            options=geneset_options, default=next(iter(geneset_options), None),
            options_callback=lambda x, y: get_geneset_options(x, y),
            master_cid=[self.children["select_groupby"].cid, self.children["select_reference"].cid],
        )

        # gsea_groupby_options = self.dataset.get_rank_genes_groups()

        # def get_gsea_reference_options():
        #     gsea_groupby_options = self.dataset.get_rank_genes_groups()
        #     if len(gsea_groupby_options) > 0:
        #         gsea_refs = sorted(list(self.dataset.adata.uns["de"][gsea_groupby_options[0]].keys()))
        #     else:
        #         gsea_refs = []
            
        #     return gsea_refs

        # gsea_refs = get_gsea_reference_options()

        # self.children.update(
        #     select_gsea_groupby=DropDown(
        #         cid=CID(page_id, loc_class, "select-gsea_groupby"),
        #         options=gsea_groupby_options, default=next(iter(gsea_groupby_options), None),
        #         options_callback=lambda: self.dataset.get_rank_genes_groups(),
        #         update_store_id="de-store"
        #     ),
        #     select_gsea_reference=DropDown(
        #         cid=CID(page_id, loc_class, "select-gsea_reference"),
        #         options=gsea_refs, default=next(iter(gsea_refs), None),
        #         options_callback=lambda: get_gsea_reference_options(),
        #         update_store_id="de-store"
        #     )
        # )

        # gene_sets = gseapy.get_library_name()
        # default_gene_set = "KEGG_2021_Human" if "KEGG_2021_Human" in gene_sets else next(iter(gene_sets), None)

        # if len(groupby_options) > 0:
        #     ref_options = get_reference_options(groupby_options[0])
        # else:
        #     ref_options = []

        # self.children["select_library"] = DropDown(
        #     cid=CID(page_id, loc_class, "select-library"),
        #     options=gene_sets, default=default_gene_set,
        # )

        self.children.update(
            term_card=TermCard(
                self.page_id, self.loc_class, self.dataset,
                select_groupby_cid=self.children["select_groupby"].cid,
                select_reference_cid=self.children["select_reference"].cid,
                select_library_cid=self.children["select_library"].cid
            )
        )
        self.actions.update(
            # apply_gsea=ApplyGSEA(
            #     CID(page_id, loc_class, "apply_gsea"), dataset,
            #     select_gsea_groupby_cid=self.children["select_gsea_groupby"].cid,
            #     select_gsea_reference_cid=self.children["select_gsea_reference"].cid,
            #     select_geneset_cid=self.children["select_library"].cid,
            # ),
            plot=Plot(
                self.cid, dataset,
                self.children["select_groupby"].cid,
                self.children["select_reference"].cid,
                self.children["select_library"].cid,
            ),
            # select_term=TermComponents.SelectTerm(CID(page_id, loc_class, "select_term"), dataset),
        )


    def create_layout(self):
        select_tab = components.FigureHeaderTab(self.page_id, self.loc_class, tab_label="Select", content=[
            html.Div([
                html.Label("Group By"),
                self.children["select_groupby"].create_layout(),
                self.children["select_groupby"].get_stores()
            ], className="param-row-stacked"),
            # Projection Hue celect
            html.Div([
                html.Label("Reference"),
                self.children["select_reference"].create_layout(),
                self.children["select_reference"].get_stores()
            ], className="param-row-stacked"),
            html.Div([
                html.Label("Library"),
                self.children["select_library"].create_layout(),
                self.children["select_library"].get_stores()
            ], className="param-row-stacked")
        ])

        term_tab = components.FigureHeaderTab(
            self.page_id, self.loc_class, tab_label="Term", id=f"{self.page_id}-{self.loc_class}-termcard",
            content=self.children["term_card"].create_layout()
        )

        header = components.FigureHeader(self.page_id, self.loc_class, tabs=[select_tab, term_tab])

        figure = html.Div(
            children=[
                html.Div(children=header.create_layout(), className="fig-header"),

                html.Div(
                    [
                        dcc.Loading(type="circle", children=[
                            html.Div(
                                dcc.Graph(id=self.cid.to_dict(), className=f"{self.loc_class}-plot")
                            )],
                        )
                    ],
                    className=f"{self.loc_class}-body",
                ),
            ],
            className=f"{self.loc_class}",
        )
        return figure
    
    def get_sidebar_params(self) -> list:
        divs = [
            # html.Div(
            #     children=[
            #         html.Label("GSEA Group By"),
            #         self.children["select_gsea_groupby"].create_layout(),
            #         self.children["select_gsea_groupby"].get_stores()
            #     ],
            #     className="param-row-stacked",
            # ),
            # html.Div(
            #     children=[
            #         html.Label("GSEA Reference"),
            #         self.children["select_gsea_reference"].create_layout(),
            #         self.children["select_gsea_reference"].get_stores()
            #     ],
            #     className="param-row-stacked",
            # ),
            # html.Div(
            #     children=[
            #         html.Label("Library"),
            #         self.children["select_library"].create_layout(),
            #         self.children["select_library"].get_stores()
            #     ],
            #     className="param-row-stacked",
            # )
        ]


        return divs