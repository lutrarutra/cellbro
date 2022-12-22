import dash
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html

from cellbro.util.DashAction import DashAction
from cellbro.util.DashPage import DashPage
import cellbro.util.Components as Components
import cellbro.plots.Heatmap as Heatmap

import scout

gsea_volcano_layout = dict(
    paper_bgcolor="white",
    plot_bgcolor="white",
    # xaxis=dict(showgrid=False, zeroline=False,visible=True, showticklabels=True),
    # yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=5, b=5, l=5, r=5),
)

class ListAvailableLRefs(DashAction):
    def setup_callbacks(self, app):
        output = [
            Output(component_id="gsea-groupby", component_property="options"),
            Output(component_id="gsea-reference", component_property="options"),
        ]
        inputs = {
            "groupby": Input(component_id="gsea-groupby", component_property="value"),
            "store": Input(component_id="groupby-store", component_property="data"),
        }

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(groupby, store):
            if groupby is None and store is None:
                raise PreventUpdate

            rgg = self.dataset.get_rank_genes_groups()
            return [
                rgg, sorted(list(self.dataset.adata.uns[f"rank_genes_{rgg[0]}"].keys()))
            ]
        

class ApplyGSEA(DashAction):
    def setup_callbacks(self, app):
        output = [
            Output(component_id="gsea-volcano-plot", component_property="figure"),
        ]
        inputs = {
            "submit": Input(component_id="gsea-submit", component_property="n_clicks"),
        }
        state = {
            "groupby": State(component_id="gsea-groupby", component_property="value"),
            "reference": State(component_id="gsea-reference", component_property="value"),
        }

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, groupby, reference):
            if submit is None:
                raise PreventUpdate

            gene_score_df = self.dataset.adata.uns[f"rank_genes_{groupby}"][reference]
            res = scout.tl.GSEA(gene_score_df, score_of_interest="gene_score")
            self.dataset.adata.uns[f"gsea_{groupby}_{reference}"] = res
            
            return [scout.ply.gsea_volcano(res, layout=gsea_volcano_layout)]


class GSEAPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.gsea", "GSEA", "gsea", order)
        self.dataset = dataset
        self.heatmap = Heatmap.Heatmap(dataset, self.id)
        self.actions.update(
            apply_de=ApplyGSEA(dataset=self.dataset, page_id_prefix=self.id),
            list_available_refs=ListAvailableLRefs(dataset=self.dataset, page_id_prefix=self.id),
        )

    def create_layout(self) -> list:
        top_sidebar = Components.create_sidebar(
            id="gsea-top-sidebar", btn_id="gsea-submit",
            title="Gene Set Enrichment Settings",
            params_children=self._params_layout(),
            class_name="top-sidebar"
        )

        main = html.Div(
            children=[
                html.Div(
                    children=[
                        # Components.create_gene_card(None, self.dataset),
                    ],
                    id="gsea-volcano-info",
                    className="main-select top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            id="loading-gsea-volcano",
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id="gsea-volcano-plot", className="main-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    id="gsea-volcano-figure",
                    className="main-figure",
                ),
            ],
            className="main",
        )

        bottom_left_sidebar, bottom_figure = self.heatmap.create_layout()

        bot_sidebar = Components.create_sidebar(
            id="gsea-bot-sidebar", btn_id=None,
            title="Gene Set Enrichment Settings",
            params_children=[],
            class_name="bot-sidebar"
        )

        layout = [
            html.Div(
                id="top", className="top", children=[top_sidebar, main]
            ),
            html.Div(
                id="bottom",
                className="bottom",
                children=[
                    # bottom_sidebar, bottom_figure
                    bot_sidebar
                ],
            ),
        ]
        return layout

    def _params_layout(self):
        rank_genes_groups = self.dataset.get_rank_genes_groups()
        if len(rank_genes_groups) > 0:
            gsea_refs = sorted(list(self.dataset.adata.uns[f"rank_genes_{rank_genes_groups[0]}"].keys()))
        else:
            gsea_refs = []

        divs = [
            html.Div(
                children=[
                    html.Label("GSEA GroupBy"),
                    html.Div(
                         [
                            dcc.Dropdown(
                                options=rank_genes_groups,
                                value=next(iter(rank_genes_groups), None),
                                id="gsea-groupby", clearable=False,
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
                                id="gsea-reference", clearable=False,
                            ),
                        ],
                        className="param-select",
                    ),
                ],
                className="param-row-stacked",
            )
        ]
        return divs
