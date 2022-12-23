import dash
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, ctx

from cellbro.util.DashAction import DashAction
from cellbro.util.DashPage import DashPage
import cellbro.util.Components as Components
import cellbro.plots.Heatmap as Heatmap
import cellbro.plots.Projection as Projection

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
            Output(component_id=f"{self.page_id_prefix}-groupby", component_property="options"),
            Output(component_id=f"{self.page_id_prefix}-groupby", component_property="value"),
            Output(component_id=f"{self.page_id_prefix}-reference", component_property="options"),
            Output(component_id=f"{self.page_id_prefix}-reference", component_property="value"),
        ]
        inputs = {
            "groupby": Input(component_id=f"{self.page_id_prefix}-groupby", component_property="value"),
        }

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(groupby):
            rank_genes_groups = self.dataset.get_rank_genes_groups()
            if len(rank_genes_groups) == 0:
                raise PreventUpdate

            if groupby is None:
                refs = sorted(list(self.dataset.adata.uns[f"rank_genes_{rank_genes_groups[0]}"].keys()))
            else:
                refs = sorted(list(self.dataset.adata.uns[f"rank_genes_{groupby}"].keys()))

            return [
                rank_genes_groups, rank_genes_groups[0],
                refs, refs[0]
            ]
        

class ApplyGSEA(DashAction):
    def apply(self, groupby, reference):
        gene_score_df = self.dataset.adata.uns[f"rank_genes_{groupby}"][reference]
        res = scout.tl.GSEA(gene_score_df, score_of_interest="gene_score")
        self.dataset.adata.uns[f"gsea_{groupby}_{reference}"] = res

        return [scout.ply.gsea_volcano(res, layout=gsea_volcano_layout)]

    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-volcano-plot", component_property="figure"),
        ]
        inputs = {
            "submit": Input(component_id=f"{self.page_id_prefix}-submit", component_property="n_clicks"),
        }
        state = {
            "groupby": State(component_id=f"{self.page_id_prefix}-groupby", component_property="value"),
            "reference": State(component_id=f"{self.page_id_prefix}-reference", component_property="value"),
        }

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, groupby, reference):
            if submit is None:
                raise PreventUpdate

            return self.apply(groupby, reference)


class PlotHeatmap(DashAction):
    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-heatmap-selected_genes", component_property="value"),
            Output(f"{self.page_id_prefix}-heatmap-cluster_cells_by", "value"),
            Output(component_id=f"{self.page_id_prefix}-heatmap-plot", component_property="figure"),
            Output(component_id=f"{self.page_id_prefix}-heatmap-plot", component_property="style"),
        ]

        # Inputs to Projection
        inputs = {
            "submit": Input(
                component_id=f"{self.page_id_prefix}-heatmap-submit", component_property="n_clicks"
            ),
            "click_data": Input(f"{self.page_id_prefix}-volcano-plot", "clickData"),
        }

        state = dict(
            selected_genes=State(f"{self.page_id_prefix}-heatmap-selected_genes", "value"),
            cluster_cells_by=State(f"{self.page_id_prefix}-heatmap-cluster_cells_by", "value"),
            gsea_groupby=State(component_id=f"{self.page_id_prefix}-groupby", component_property="value"),
            gsea_reference=State(component_id=f"{self.page_id_prefix}-reference", component_property="value"),
            categoricals=State(f"{self.page_id_prefix}-heatmap-selected_categoricals", "value")

        )
        for key in Heatmap.heatmap_params.keys():
            state[key] = State(component_id=f"{self.page_id_prefix}-heatmap-{key}", component_property="value")

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, click_data, **kwargs):
            if click_data is None:
                raise PreventUpdate
            
            selected_genes = kwargs["selected_genes"]
            cluster_cells_by = kwargs["cluster_cells_by"]

            if ctx.triggered_id == f"{self.page_id_prefix}-volcano-plot":
                cluster_cells_by = kwargs["gsea_groupby"]
                reference = kwargs["gsea_reference"]

                term = click_data["points"][0]["hovertext"]
                res = self.dataset.adata.uns[f"gsea_{cluster_cells_by}_{reference}"]
                selected_genes = res[res["Term"] == term]["lead_genes"].values[0]
                
                kwargs["selected_genes"] = selected_genes
                kwargs["cluster_cells_by"] = cluster_cells_by

            fig, style = Heatmap.Heatmap.plot(self.dataset, kwargs)

            return [selected_genes, cluster_cells_by, fig, style]

class GSEAPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.gsea", "GSEA", "gsea", order)
        self.dataset = dataset
        self.actions.update(
            apply_de=ApplyGSEA(dataset=self.dataset, page_id_prefix=self.id),
            list_available_refs=ListAvailableLRefs(dataset=self.dataset, page_id_prefix=self.id),
        )
        self.plots.update(
            projection=Projection.Projection(dataset, self.id, loc_class="secondary"),
            heatmap=Heatmap.Heatmap(dataset, self.id, loc_class="bottom")
        )

        self.plots["heatmap"].actions["plot_heatmap"] = PlotHeatmap(dataset, self.id)
        # TODO: Fix this
        self.plots["heatmap"].actions.pop("add_genes_from_list")

    def create_layout(self) -> list:
        self.sidebars["left_sidebar"] = Components.Sidebar(
            page_id_prefix=self.id, apply_btn_id=f"{self.id}-submit",
            title="Gene Set Enrichment Settings",
            params_children=self._params_layout(),
            row="top", side="left",
        )

        main = html.Div(
            children=[
                html.Div(
                    children=[
                        # Components.create_gene_card(None, self.dataset),
                    ],
                    id=f"{self.id}-volcano-info",
                    className="main-select top-parameters",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id=f"{self.id}-volcano-plot", className="main-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    id=f"{self.id}-volcano-figure",
                    className="main-figure",
                ),
            ],
            className="main",
        )

        self.sidebars["right_sidebar"] = Components.Sidebar(
            page_id_prefix=self.id, apply_btn_id=f"{self.id}-projection-submit",
            title="Projection Settings",
            params_children=self.plots["projection"].get_sidebar_params(),
            row="top", side="right",
        )
        secondary_figure = self.plots["projection"].create_layout()
        
        self.sidebars["bot_sidebar"] = Components.Sidebar(
            page_id_prefix=self.id, row="bot", side="left",
            title="Heatmap Settings",
            params_children=self.plots["heatmap"].get_sidebar_params(),
            apply_btn_id=f"{self.id}-heatmap-submit", btn_text="Plot"
        )
        
        bot_figure = self.plots["heatmap"].create_layout()

        layout = [
            html.Div(
                id="top", className="top", children=[
                    self.sidebars["left_sidebar"].create_layout(), main,
                    self.sidebars["right_sidebar"].create_layout(), secondary_figure
                ]
            ),
            html.Div(
                id="bottom",
                className="bottom",
                children=[
                    # bottom_sidebar, bottom_figure
                    self.sidebars["bot_sidebar"].create_layout(), bot_figure
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
                    html.Label("GSEA Group By"),
                    html.Div(
                         [
                            dcc.Dropdown(
                                options=rank_genes_groups,
                                value=next(iter(rank_genes_groups), None),
                                id=f"{self.id}-groupby", clearable=False,
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
                                id=f"{self.id}-reference", clearable=False,
                            ),
                        ],
                        className="param-select",
                    ),
                ],
                className="param-row-stacked",
            )
        ]
        return divs
