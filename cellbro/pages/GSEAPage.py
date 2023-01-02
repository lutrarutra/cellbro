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

import scout

class ListAvailableLRefs(DashAction):
    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id}-gsea-groupby", component_property="options"),
            Output(component_id=f"{self.page_id}-gsea-groupby", component_property="value"),
            Output(component_id=f"{self.page_id}-gsea-reference", component_property="options"),
            Output(component_id=f"{self.page_id}-gsea-reference", component_property="value"),
        ]
        inputs = {
            "groupby": Input(component_id=f"{self.page_id}-gsea-groupby", component_property="value"),
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


class PlotHeatmap(DashAction):
    RType = namedtuple("RType", ["figure", "style"])

    def setup_callbacks(self, app):
        output = [
            # Output(f"{self.page_id}-{self.loc_class}-param-genes", "value"),
            # Output(f"{self.page_id}-{self.loc_class}-param-cluster_cells_by", "value"),
            Output(f"{self.page_id}-heatmap-plot", "figure"),
            Output(f"{self.page_id}-heatmap-plot", "style"),
        ]

        # Inputs to Projection
        inputs = dict(
            submit=Input(f"{self.page_id}-{self.loc_class}-sidebar-apply_btn", "n_clicks"),
            click_data=Input(f"{self.page_id}-main-plot", "clickData"),
            gsea_groupby = Input(f"{self.page_id}-main-groupby", "value"),
            gsea_reference = Input(f"{self.page_id}-main-reference", "value"),
        )

        state = dict(
            selected_genes=State(f"{self.page_id}-{self.loc_class}-param-genes", "value"),
            cluster_cells_by=State(f"{self.page_id}-{self.loc_class}-param-cluster_cells_by", "value"),
            categoricals=State(f"{self.page_id}-{self.loc_class}-param-categoricals", "value")
        )
        for key in heatmap_params.keys():
            state[key] = State(f"{self.page_id}-{self.loc_class}-param-{key}", "value")

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, click_data, **kwargs):
            if click_data is None:
                raise PreventUpdate
            
            selected_genes = kwargs["selected_genes"]
            cluster_cells_by = kwargs["cluster_cells_by"]

            if ctx.triggered_id == f"{self.page_id}-main-plot":
                cluster_cells_by = kwargs["gsea_groupby"]
                reference = kwargs["gsea_reference"]

                term = click_data["points"][0]["hovertext"]
                res = self.dataset.adata.uns[f"gsea"][cluster_cells_by][reference]
                selected_genes = res[res["Term"] == term]["lead_genes"].values[0]
                
                kwargs["selected_genes"] = selected_genes
                kwargs["cluster_cells_by"] = cluster_cells_by

            fig, style = Heatmap.Heatmap.plot(self.dataset, kwargs)

            return self.RType(
                # selected_genes=selected_genes,
                # cluster_cells_by=cluster_cells_by,
                figure=fig,
                style=style
            )

class PlotProjection(DashAction):
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
        outputs = [
            Output(f"{self.page_id}-projection-plot", "figure"),
            Output(f"{self.page_id}-{self.loc_class}-select-projection_type", "options"),
            Output(f"{self.page_id}-{self.loc_class}-select-projection_type", "value"),
        ]

        inputs = dict(
            submit=Input(f"{self.page_id}-main-sidebar-apply_btn", "n_clicks"),
            color=Input(f"{self.page_id}-{self.loc_class}-select-color", "value"),
            obsm_layer=Input(f"{self.page_id}-{self.loc_class}-select-projection_type", "value"),
            continuous_cmap=Input(f"{self.page_id}-{self.loc_class}-select-continuous_cmap", "value"),
            discrete_cmap=Input(f"{self.page_id}-{self.loc_class}-select-discrete_cmap", "value"),
            click_data=Input(f"{self.page_id}-main-plot", "clickData"),
        )

        state = dict(
            projection_type=State(f"{self.page_id}-{self.loc_class}-type-select", "value"),
            gsea_groupby=State(f"{self.page_id}-{self.loc_class}-groupby", "value"),
            gsea_reference=State(f"{self.page_id}-{self.loc_class}-reference", "value"),
        )
        for key in prj.UMAP._params.keys():
            state[f"umap_{key}"] = State(f"{self.page_id}-projection-umap-{key}", "value")

        for key in prj.SCVI_UMAP._params.keys():
            state[f"scvi_umap_{key}"] = State(f"{self.page_id}-projection-scvi_umap-{key}", "value")

        for key in prj.TSNE._params.keys():
            state[f"tsne_{key}"] = State(f"{self.page_id}-projection-tsne-{key}", "value")

        for key in prj.Trimap._params.keys():
            state[f"trimap_{key}"] = State(f"{self.page_id}-projection-trimap-{key}", "value")

        @app.dash_app.callback(output=outputs, inputs=inputs, state=state)
        def _(submit, color, obsm_layer, projection_type, continuous_cmap, discrete_cmap, click_data, **kwargs):
            if ctx.triggered_id == f"{self.page_id}-projection-submit":
                if projection_submit is not None:
                    obsm_layer = self.apply(projection_type, params=kwargs)

            if click_data is not None:
                groupby = kwargs["gsea_groupby"]
                reference = kwargs["gsea_reference"]
                term = click_data["points"][0]["hovertext"]
                res = self.dataset.adata.uns[f"gsea"][groupby][reference]
                selected_genes = res[res["Term"] == term]["lead_genes"].values[0]
                fig = self.plot(
                    color=selected_genes, obsm_layer=obsm_layer, hue_aggregate=None,
                    continuous_cmap=continuous_cmap, discrete_cmap=discrete_cmap
                )
            else:
                fig = self.plot(
                    color=color, obsm_layer=obsm_layer,
                    continuous_cmap=continuous_cmap, discrete_cmap=discrete_cmap
                )

            return (fig, list(self.dataset.adata.obsm.keys()), obsm_layer)

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
            heatmap=Heatmap(dataset, self.page_id, loc_class="bottom")
        )

        self.components["heatmap"].actions["plot_heatmap"] = PlotHeatmap(CID(self.page_id, "bottom", "plot_heatmap"), self.dataset)
        # # TODO: Fix this
        # self.components["heatmap"].actions.pop("edit_genelist")
        self.components["projection"].actions["plot_projection"] = PlotProjection(CID(self.page_id, "secondary", "plot_projection"), dataset)

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
        
        self.components["bot_sidebar"] = Sidebar(
            page_id=self.page_id, loc_class="bottom",
            title="Heatmap Settings", create_btn=True, btn_text="Plot",
            params_children=self.components["heatmap"].get_sidebar_params(),
        )
        
        bot_figure = self.components["heatmap"].create_layout()

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
                    self.components["bot_sidebar"].create_layout(), bot_figure
                ],
            ),
        ]
        return layout
