from functools import namedtuple

import dash
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, ctx

from ..util.DashAction import DashAction
from ..components.DashPage import DashPage
from ..plots.Heatmap import Heatmap, heatmap_params
from ..plots import projection as prj
from ..plots.GSEA.GSEAVolcano import GSEAVolcano
from ..components.CID import CID
from ..components.Sidebar import Sidebar
from ..components import components

import scout

# class ListAvailableLRefs(DashAction):
#     def setup_callbacks(self, app):
#         output = [
#             Output(component_id=f"{self.page_id}-gsea-groupby", component_property="options"),
#             Output(component_id=f"{self.page_id}-gsea-groupby", component_property="value"),
#             Output(component_id=f"{self.page_id}-gsea-reference", component_property="options"),
#             Output(component_id=f"{self.page_id}-gsea-reference", component_property="value"),
#         ]
#         inputs = {
#             "groupby": Input(component_id=f"{self.page_id}-gsea-groupby", component_property="value"),
#         }

#         @app.dash_app.callback(output=output, inputs=inputs)
#         def _(groupby):
#             rank_genes_groups = self.dataset.get_rank_genes_groups()
#             if len(rank_genes_groups) == 0:
#                 raise PreventUpdate

#             if groupby is None:
#                 refs = sorted(list(self.dataset.adata.uns[f"rank_genes_{rank_genes_groups[0]}"].keys()))
#             else:
#                 refs = sorted(list(self.dataset.adata.uns[f"rank_genes_{groupby}"].keys()))

#             return [
#                 rank_genes_groups, rank_genes_groups[0],
#                 refs, refs[0]
#             ]


class PlotHeatmap(DashAction):
    RType = namedtuple("RType", ["figure", "style"])

    def setup_callbacks(self, app):
        output = [
            Output(f"{self.page_id}-heatmap-plot", "figure"),
            Output(f"{self.page_id}-heatmap-plot", "style"),
        ]

        # Inputs to Projection
        inputs = dict(
            submit=Input(f"{self.page_id}-{self.loc_class}-sidebar-apply_btn", "n_clicks"),
            click_data=Input(f"{self.page_id}-main-plot", "clickData"),
            gsea_groupby = Input(f"{self.page_id}-main-groupby", "value"),
            gsea_reference = Input(f"{self.page_id}-main-reference", "value"),
            selected_genes=State(f"{self.page_id}-{self.loc_class}-param-genes", "value"),
            cluster_cells_by=State(f"{self.page_id}-{self.loc_class}-param-cluster_cells_by", "value"),
            categoricals=State(f"{self.page_id}-{self.loc_class}-param-categoricals", "value")
        )

        for key in heatmap_params.keys():
            inputs[key] = State(f"{self.page_id}-{self.loc_class}-param-{key}", "value")

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(submit, click_data, gsea_groupby, gsea_reference, selected_genes, cluster_cells_by, categoricals, **kwargs):
            if click_data is None:
                raise PreventUpdate

            if ctx.triggered_id == f"{self.page_id}-main-plot":
                cluster_cells_by = gsea_groupby
                reference = gsea_reference

                term = click_data["points"][0]["hovertext"]
                res = self.dataset.adata.uns[f"gsea"][cluster_cells_by][reference]
                selected_genes = res[res["Term"] == term]["lead_genes"].values[0]
                
                selected_genes = selected_genes
                cluster_cells_by = cluster_cells_by

            fig, style = Heatmap.Heatmap.plot(self.dataset, kwargs)

            return self.RType(
                # selected_genes=selected_genes,
                # cluster_cells_by=cluster_cells_by,
                figure=fig,
                style=style
            )

class PlotProjection(DashAction):
    def __init__(
        self, parent_cid: CID, dataset,
        obsm_layer_cid: CID,
        continuous_cmap_cid: CID,
        discrete_cmap_cid: CID,
        gsea_volcano_cid: CID,
        select_groupby_cid: CID,
        select_reference_cid: CID,
    ):
        super().__init__(parent_cid, dataset)
        self.obsm_layer_cid = obsm_layer_cid
        self.continuous_cmap_cid = continuous_cmap_cid
        self.discrete_cmap_cid = discrete_cmap_cid
        self.gsea_volcano_cid = gsea_volcano_cid
        self.select_groupby_cid = select_groupby_cid
        self.select_reference_cid = select_reference_cid

    def plot(self, color, obsm_layer, continuous_cmap, discrete_cmap, **kwargs):
        fig = scout.ply.projection(
            self.dataset.adata, obsm_layer=obsm_layer, hue=color,
            layout=prj.prj_tools.default_layout, 
            continuous_cmap=continuous_cmap, discrete_cmap=discrete_cmap,
            **kwargs
        )
        return fig

    def apply(self, projection_type, params):
        if projection_type == "UMAP":
            projection = prj.UMAP(self.dataset, prj.UMAP.parse_params(params))
        elif projection_type == "t-SNE":
            projection = prj.TSNE(self.dataset, prj.TSNE.parse_params(params))
        elif projection_type == "Trimap":
            projection = prj.Trimap(self.dataset, prj.Trimap.parse_params(params))

        return projection.apply()

    def setup_callbacks(self, app):
        outputs = Output(self.parent_cid.to_dict(), "figure"),

        inputs = dict(
            obsm_layer=Input(self.obsm_layer_cid.to_dict(), "value"),
            continuous_cmap=Input(self.continuous_cmap_cid.to_dict(), "value"),
            discrete_cmap=Input(self.discrete_cmap_cid.to_dict(), "value"),
            click_data=Input(self.gsea_volcano_cid.to_dict(), "clickData"),
            groupby=State(self.select_groupby_cid.to_dict(), "value"),
            reference=State(self.select_reference_cid.to_dict(), "value"),
        )
        print("YAHOO")
        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(obsm_layer, continuous_cmap, discrete_cmap, click_data, groupby, reference):
            print(ctx.triggered_id)
            if click_data is None:
                raise PreventUpdate

            term = click_data["points"][0]["hovertext"]
            res = self.dataset.adata.uns[f"gsea"][groupby][reference]
            selected_genes = res[res["Term"] == term]["lead_genes"].values[0]
            fig = self.plot(
                color=selected_genes, obsm_layer=obsm_layer, hue_aggregate=None,
                continuous_cmap=continuous_cmap, discrete_cmap=discrete_cmap
            )
            return fig

class GSEAPage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.gsea", "GSEA", "gsea", order)
        self.dataset = dataset
        self.actions.update(
            # list_available_refs=ListAvailableLRefs(CID(self.page_id, "static", "list_available_refs"), self.dataset),
        )
        self.components.update(
            main=GSEAVolcano(dataset, self.page_id, loc_class="main"),
            projection=prj.Projection(dataset, self.page_id, loc_class="secondary"),
            # heatmap=Heatmap(dataset, self.page_id, loc_class="bottom")
        )

        # self.components["heatmap"].actions["plot_heatmap"] = PlotHeatmap(CID(self.page_id, "bottom", "plot_heatmap"), self.dataset)
        # # TODO: Fix this
        # self.components["heatmap"].actions.pop("edit_genelist")

        self.components["projection"].actions["plot_projection"] = PlotProjection(
            self.components["projection"].cid, dataset,
            obsm_layer_cid=self.components["projection"].children["select_projection_type"].cid,
            continuous_cmap_cid=self.components["projection"].children["select_continuous_cmap"].cid,
            discrete_cmap_cid=self.components["projection"].children["select_discrete_cmap"].cid,
            gsea_volcano_cid=self.components["main"].cid,
            select_groupby_cid=self.components["main"].children["select_groupby"].cid,
            select_reference_cid=self.components["main"].children["select_reference"].cid,
        )

    def create_layout(self) -> list:
        self.components["left_sidebar"] = Sidebar(
            page_id=self.page_id, loc_class="main",
            create_btn=True, title="Gene Set Enrichment Settings",
            params_children=self.components["main"].get_sidebar_params()
        )

        self.components["right_sidebar"] = Sidebar(
            page_id=self.page_id, loc_class="secondary",
            create_btn=True, title="Projection Settings",
            params_children=self.components["projection"].get_sidebar_params(),
        )
        secondary_figure = self.components["projection"].create_layout()
        
        # self.components["bot_sidebar"] = Sidebar(
        #     page_id=self.page_id, loc_class="bottom",
        #     title="Heatmap Settings", create_btn=True, btn_text="Plot",
        #     params_children=self.components["heatmap"].get_sidebar_params(),
        # )
        
        #bot_figure = self.components["heatmap"].create_layout()

        layout = [
            html.Div(
                id="top", className="top", children=[
                    self.components["left_sidebar"].create_layout(), self.components["main"].create_layout(),
                    self.components["right_sidebar"].create_layout(), secondary_figure
                ]
            ),
            html.Div(
                id="bottom",
                className="bottom",
                children=[
                    # bottom_sidebar, bottom_figure
                    #self.components["bot_sidebar"].create_layout(), bot_figure
                ],
            ),
        ]
        return layout
