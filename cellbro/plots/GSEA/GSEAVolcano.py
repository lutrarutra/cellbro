import dash
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, ctx

import gseapy

from ...components.DashFigure import DashFigure
from ...util.DashAction import DashAction
from ...components import components
from ...components import TermComponents

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
            Output(component_id=f"{self.page_id_prefix}-{self.loc_class}-plot", component_property="figure"),
            Output(component_id=f"{self.page_id_prefix}-{self.loc_class}-groupby", component_property="options"),
            Output(component_id=f"{self.page_id_prefix}-{self.loc_class}-groupby", component_property="value"),
            Output(component_id=f"{self.page_id_prefix}-{self.loc_class}-reference", component_property="options"),
            Output(component_id=f"{self.page_id_prefix}-{self.loc_class}-reference", component_property="value")
        ]
        inputs = dict(
            submit=Input(component_id=f"{self.page_id_prefix}-submit", component_property="n_clicks"),
            groupby=State(component_id=f"{self.page_id_prefix}-{self.loc_class}-groupby", component_property="value"),
            reference=State(component_id=f"{self.page_id_prefix}-{self.loc_class}-reference", component_property="value"),
        )

        state = dict(
            de_groupby=State(component_id=f"{self.page_id_prefix}-gsea-groupby", component_property="value"),
            de_reference=State(component_id=f"{self.page_id_prefix}-gsea-reference", component_property="value"),
            de_gene_set=State(component_id=f"{self.page_id_prefix}-gsea-gene_set", component_property="value"),
        )

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, groupby, reference, de_groupby, de_reference, de_gene_set):
            if ctx.triggered_id == f"{self.page_id_prefix}-submit":
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



class GSEAVolcano(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.actions.update(
            apply_gsea=ApplyGSEA(dataset, page_id_prefix, loc_class),
            select_term=TermComponents.SelectTerm(dataset, page_id_prefix, loc_class),
        )

    def create_layout(self):
        select_tab = components.FigureHeaderTab(self.page_id_prefix, tab_label="Select", children=[
            html.Div([
                html.Label("Group By"),
                dcc.Dropdown(
                    options=[], value=None,
                    id=f"{self.page_id_prefix}-{self.loc_class}-groupby", clearable=False,
                ),
            ], className="param-row-stacked"),
            # Projection Hue celect
            html.Div([
                html.Label("Color"),
                dcc.Dropdown(
                    options=[], value=None,
                    id=f"{self.page_id_prefix}-{self.loc_class}-reference",
                    clearable=False,
                )
            ], className="param-row-stacked")
        ])

        term_tab = components.FigureHeaderTab(
            self.page_id_prefix, tab_label="Term", id=f"{self.page_id_prefix}-{self.loc_class}-termcard",
            children=[
                TermComponents.create_term_card(None, None, None)
            ]
        )

        header = components.FigureHeader(self.page_id_prefix, tabs=[select_tab, term_tab])

        figure = html.Div(
            children=[
                html.Div(children=header.create_layout(), className="fig-header"),

                html.Div(
                    [
                        dcc.Loading(type="circle", children=[
                            html.Div(
                                dcc.Graph(
                                    id=f"{self.page_id_prefix}-{self.loc_class}-plot", className=f"{self.loc_class}-plot"
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
                                id=f"{self.page_id_prefix}-gsea-groupby", clearable=False,
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
                                id=f"{self.page_id_prefix}-gsea-reference", clearable=False,
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
                                id=f"{self.page_id_prefix}-gsea-gene_set", clearable=False,
                            ),
                        ],
                        className="param-select",
                    ),
                ],
                className="param-row-stacked",
            )
        ]


        return divs