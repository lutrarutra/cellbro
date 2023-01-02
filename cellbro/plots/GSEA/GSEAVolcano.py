import dash
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, ctx

import gseapy

from ...components.DashPlot import DashPlot
from ...util.DashAction import DashAction
from ...components import components
from ...components import TermComponents
from ...components.CID import CID

import scout

default_layout = dict(
    paper_bgcolor="white",
    plot_bgcolor="white",
    margin=dict(t=5, b=5, l=5, r=5),
)

class ApplyGSEA(DashAction):
    def apply(self, groupby, reference, gene_set):
        gene_score_df = self.dataset.adata.uns[f"rank_genes_{groupby}"][reference]
        res = scout.tl.GSEA(gene_score_df, score_of_interest="gene_score", gene_set=gene_set)
        self.dataset.add_gsea_result(res, groupby, reference)

    def plot(self, groupby, reference):
        res = self.dataset.adata.uns[f"gsea"][groupby][reference]
        return scout.ply.gsea_volcano(res, layout=default_layout)

    def setup_callbacks(self, app):
        output = [
            Output(f"{self.page_id}-{self.loc_class}-plot", "figure"),
            Output(f"{self.page_id}-{self.loc_class}-groupby", "options"),
            Output(f"{self.page_id}-{self.loc_class}-groupby", "value"),
            Output(f"{self.page_id}-{self.loc_class}-reference", "options"),
            Output(f"{self.page_id}-{self.loc_class}-reference", "value")
        ]
        inputs = dict(
            submit=Input(f"{self.page_id}-{self.loc_class}-sidebar-apply_btn", "n_clicks"),
            groupby=State(f"{self.page_id}-{self.loc_class}-groupby", "value"),
            reference=State(f"{self.page_id}-{self.loc_class}-reference", "value"),
        )

        state = dict(
            de_groupby=State(f"{self.page_id}-gsea-groupby", "value"),
            de_reference=State(f"{self.page_id}-gsea-reference", "value"),
            de_gene_set=State(f"{self.page_id}-gsea-gene_set", "value"),
        )

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, groupby, reference, de_groupby, de_reference, de_gene_set):
            if ctx.triggered_id == f"{self.page_id}-{self.loc_class}-sidebar-apply_btn":
                if submit is not None:
                    self.apply(de_groupby, de_reference, de_gene_set)
                    groupby = de_groupby
                    reference = de_reference
                else:
                    raise PreventUpdate

            if "gsea" not in self.dataset.adata.uns:
                return [dash.no_update, [], None, [], None]

            groupby_options = list(self.dataset.adata.uns["gsea"].keys())
            if groupby is None:
                if len(groupby_options) == 0:
                    return [dash.no_update, [], None, [], None]

                groupby = groupby_options[0]

            reference_options = list(self.dataset.adata.uns[f"gsea"][groupby].keys())

            if reference is None:
                reference = reference_options[0]

            fig = self.plot(groupby, reference)

            return [fig, groupby_options, groupby, reference_options, reference]



class GSEAVolcano(DashPlot):
    def __init__(self, dataset, page_id, loc_class):
        super().__init__(dataset, page_id, loc_class)
        self.actions.update(
            apply_gsea=ApplyGSEA(CID(page_id, loc_class, "apply_gsea"), dataset),
            select_term=TermComponents.SelectTerm(CID(page_id, loc_class, "select_term"), dataset),
        )

    def create_layout(self):
        select_tab = components.FigureHeaderTab(self.page_id, self.loc_class, tab_label="Select", children=[
            html.Div([
                html.Label("Group By"),
                dcc.Dropdown(
                    options=[], value=None,
                    id=f"{self.page_id}-{self.loc_class}-groupby", clearable=False,
                ),
            ], className="param-row-stacked"),
            # Projection Hue celect
            html.Div([
                html.Label("Color"),
                dcc.Dropdown(
                    options=[], value=None,
                    id=f"{self.page_id}-{self.loc_class}-reference",
                    clearable=False,
                )
            ], className="param-row-stacked")
        ])

        term_tab = components.FigureHeaderTab(
            self.page_id, self.loc_class, tab_label="Term", id=f"{self.page_id}-{self.loc_class}-termcard",
            children=[
                TermComponents.create_term_card(None, None, None)
            ]
        )

        header = components.FigureHeader(self.page_id, self.loc_class, tabs=[select_tab, term_tab])

        figure = html.Div(
            children=[
                html.Div(children=header.create_layout(), className="fig-header"),

                html.Div(
                    [
                        dcc.Loading(type="circle", children=[
                            html.Div(
                                dcc.Graph(
                                    id=f"{self.page_id}-{self.loc_class}-plot", className=f"{self.loc_class}-plot"
                                )
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
        rank_genes_groups = self.dataset.get_rank_genes_groups()
        gene_sets = gseapy.get_library_name()
        default_gene_set = "KEGG_2021_Human" if "KEGG_2021_Human" in gene_sets else next(iter(gene_sets), None)

        if len(rank_genes_groups) > 0:
            gsea_refs = sorted(list(self.dataset.adata.uns[f"rank_genes_{rank_genes_groups[0]}"].keys()))
        else:
            gsea_refs = []

        divs = [
            html.Div(
                children=[
                    html.Label("GSEA Group By"),
                    html.Div(
                         [
                            dcc.Dropdown(
                                options=rank_genes_groups,
                                value=next(iter(rank_genes_groups), None),
                                id=f"{self.page_id}-gsea-groupby", clearable=False,
                            ),
                         ],
                        className="param-select",
                    ),
                ],
                className="param-row-stacked",
            ),
            html.Div(
                children=[
                    html.Label("GSEA Reference"),
                    html.Div(
                        [
                            dcc.Dropdown(
                                options=gsea_refs,
                                value=next(iter(gsea_refs), None),
                                id=f"{self.page_id}-gsea-reference", clearable=False,
                            ),
                        ],
                        className="param-select",
                    ),
                ],
                className="param-row-stacked",
            ),
            html.Div(
                children=[
                    html.Label("Library"),
                    html.Div(
                        [
                            dcc.Dropdown(
                                options=gene_sets,
                                value=default_gene_set,
                                id=f"{self.page_id}-gsea-gene_set", clearable=False,
                            ),
                        ],
                        className="param-select",
                    ),
                ],
                className="param-row-stacked",
            )
        ]


        return divs