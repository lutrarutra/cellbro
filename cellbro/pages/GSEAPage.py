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
from ..components.DropDown import DropDown

import scout

class PlotHeatmap(DashAction):
    RType = namedtuple("RType", ["figure", "style", "selected_genes", "cluster_by", "categoricals"])
    def __init__(
        self, parent_cid: CID, dataset,
        select_colormap_cid: CID,
        gsea_volcano_cid: CID,
        select_groupby_cid: CID,
        select_reference_cid: CID,
        select_layer_cid: CID,
        select_cluster_by_cid: CID,
        select_categoricals_cid: CID,
        select_genes_cid: CID,
    ):
        super().__init__(parent_cid, dataset)
        self.select_colormap_cid = select_colormap_cid
        self.gsea_volcano_cid = gsea_volcano_cid
        self.select_groupby_cid = select_groupby_cid
        self.select_reference_cid = select_reference_cid
        self.select_layer_cid = select_layer_cid
        self.select_cluster_by_cid = select_cluster_by_cid
        self.select_categoricals_cid = select_categoricals_cid
        self.select_genes_cid = select_genes_cid

    def plot(self, selected_genes, cluster_cells_by, categoricals, layer, colormap):
        fig = scout.ply.heatmap(
            adata=self.dataset.adata, var_names=selected_genes, categoricals=categoricals,
            layer=layer, cluster_cells_by=cluster_cells_by, layout=dict(),
            cmap=colormap,
        )

        style = {
            # "width": f"{int(z.shape[0]/2)+100}px",
            "height": fig.layout["height"],
        }
        return fig, style

    def setup_callbacks(self, app):
        output = [
            Output(self.parent_cid.to_dict(), "figure"),
            Output(self.parent_cid.to_dict(), "style"),
            Output(self.select_genes_cid.to_dict(), "value"),
            Output(self.select_cluster_by_cid.to_dict(), "value"),
            Output(self.select_categoricals_cid.to_dict(), "value"),
        ]

        inputs = dict(
            submit=Input(f"{self.page_id}-{self.loc_class}-sidebar-apply_btn", "n_clicks"),
            click_data=Input(self.gsea_volcano_cid.to_dict(), "clickData"),
            gsea_groupby=Input(self.select_groupby_cid.to_dict(), "value"),
            gsea_reference=Input(self.select_reference_cid.to_dict(), "value"),
            selected_genes=State(self.select_genes_cid.to_dict(), "value"),
            cluster_cells_by=State(self.select_cluster_by_cid.to_dict(), "value"),
            categoricals=State(self.select_categoricals_cid.to_dict(), "value"),
            colormap=State(self.select_colormap_cid.to_dict(), "value"),
            layer=State(self.select_layer_cid.to_dict(), "value"),
        )

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(
            submit, click_data, gsea_groupby, gsea_reference, selected_genes,
            cluster_cells_by, categoricals, colormap, layer
            ):
            if click_data is None:
                raise PreventUpdate

            if ctx.triggered_id == self.gsea_volcano_cid.to_dict():
                cluster_cells_by = gsea_groupby
                reference = gsea_reference

                term = click_data["points"][0]["hovertext"]
                res = self.dataset.adata.uns[f"gsea"][cluster_cells_by][reference]
                selected_genes = res[res["Term"] == term]["lead_genes"].values[0]
                
                selected_genes = selected_genes
                cluster_cells_by = cluster_cells_by

                if isinstance(categoricals, str):
                    categoricals = [categoricals]

                if cluster_cells_by not in categoricals:
                    categoricals.append(cluster_cells_by)

            fig, style = self.plot(selected_genes, cluster_cells_by, categoricals, layer, colormap)

            return self.RType(
                figure=fig,
                style=style,
                selected_genes=selected_genes if selected_genes is not None else dash.no_update,
                cluster_by=cluster_cells_by if cluster_cells_by is not None else dash.no_update,
                categoricals=categoricals if categoricals is not None else dash.no_update,
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
        select_aggregation_cid: CID,
    ):
        super().__init__(parent_cid, dataset)
        self.obsm_layer_cid = obsm_layer_cid
        self.continuous_cmap_cid = continuous_cmap_cid
        self.discrete_cmap_cid = discrete_cmap_cid
        self.gsea_volcano_cid = gsea_volcano_cid
        self.select_groupby_cid = select_groupby_cid
        self.select_reference_cid = select_reference_cid
        self.select_aggregation_cid = select_aggregation_cid

    def plot(self, color, obsm_layer, continuous_cmap, discrete_cmap, hue_aggregate):
        fig = scout.ply.projection(
            self.dataset.adata, obsm_layer=obsm_layer, hue=color,
            layout=prj.prj_tools.default_layout, 
            continuous_cmap=continuous_cmap, discrete_cmap=discrete_cmap,
            hue_aggregate=hue_aggregate,
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
            aggregation=Input(self.select_aggregation_cid.to_dict(), "value"),
            groupby=State(self.select_groupby_cid.to_dict(), "value"),
            reference=State(self.select_reference_cid.to_dict(), "value"),
        )
        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(obsm_layer, continuous_cmap, discrete_cmap, click_data, groupby, reference, aggregation):
            if click_data is None:
                raise PreventUpdate

            term = click_data["points"][0]["hovertext"]
            res = self.dataset.adata.uns[f"gsea"][groupby][reference]
            selected_genes = res[res["Term"] == term]["lead_genes"].values[0]

            if aggregation == "sum": aggregation = None

            fig = self.plot(
                color=selected_genes, obsm_layer=obsm_layer, hue_aggregate=aggregation,
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
            heatmap=Heatmap(dataset, self.page_id, loc_class="bottom"),
            select_projection_agg=DropDown(
                cid=CID(self.page_id, "secondary", "select-projection_agg"),
                options={"abs":"Absolute", "sum":"Relative"}, default="abs",
            ),
        )

        self.components["heatmap"].actions["plot_heatmap"] = PlotHeatmap(
            self.components["heatmap"].cid, self.dataset,
            gsea_volcano_cid=self.components["main"].cid,
            select_groupby_cid=self.components["main"].children["select_groupby"].cid,
            select_reference_cid=self.components["main"].children["select_reference"].cid,
            select_colormap_cid=self.components["heatmap"].children["select_colormap"].cid,
            select_cluster_by_cid=self.components["heatmap"].children["select_clusterby"].cid,
            select_layer_cid=self.components["heatmap"].children["select_layer"].cid,
            select_categoricals_cid=self.components["heatmap"].children["select_categoricals"].cid,
            select_genes_cid=self.components["heatmap"].children["select_genes"].cid,
        )

        self.components["projection"].actions["plot_projection"] = PlotProjection(
            self.components["projection"].cid, dataset,
            obsm_layer_cid=self.components["projection"].children["select_projection_type"].cid,
            continuous_cmap_cid=self.components["projection"].children["select_continuous_cmap"].cid,
            discrete_cmap_cid=self.components["projection"].children["select_discrete_cmap"].cid,
            gsea_volcano_cid=self.components["main"].cid,
            select_groupby_cid=self.components["main"].children["select_groupby"].cid,
            select_reference_cid=self.components["main"].children["select_reference"].cid,
            select_aggregation_cid=self.components["select_projection_agg"].cid,
        )
        self.components["projection"].children["type_header_tab"] = components.FigureHeaderTab(
            self.page_id, self.components["projection"].loc_class, tab_label="Type", content=[
                # Projection type celect
                html.Div([
                    html.Label("Projection Type"),
                    self.components["projection"].children["select_projection_type"].create_layout(),
                    self.components["projection"].children["select_projection_type"].get_stores(),
                ], className="param-row-stacked"),
                # Projection Hue celect
                html.Div([
                    html.Label("Aggregation"),
                    self.components["select_projection_agg"].create_layout(),
                    self.components["select_projection_agg"].get_stores(),
                ], className="param-row-stacked")
            ]
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
                    self.components["bot_sidebar"].create_layout(), bot_figure
                ],
            ),
        ]
        return layout
