import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

from ...components.DashPlot import DashPlot
from ...util.DashAction import DashAction
from ...components import components
from .. import PCA
from ...components.CID import CID
from ...components.DropDown import DropDown
from ...components.InputField import InputField

import scout

class PlotPCA(DashAction):
    def __init__(
        self, parent_cid: CID, dataset,
        select_hue_cid: CID,
        select_pcx_cid: CID,
        select_pcy_cid: CID,
        select_discrete_cmap_cid: CID,
        select_continuous_cmap_cid: CID,
    ):
        super().__init__(parent_cid, dataset)
        self.select_hue_cid = select_hue_cid
        self.select_pcx_cid = select_pcx_cid
        self.select_pcy_cid = select_pcy_cid
        self.select_discrete_cmap_cid = select_discrete_cmap_cid
        self.select_continuous_cmap_cid = select_continuous_cmap_cid

    def plot(self, color, pc_x, pc_y, continuous_cmap, discrete_cmap):
        fig = scout.ply.projection(
            self.dataset.adata, obsm_layer="X_pca", hue=color, components=[pc_x, pc_y],
            layout=PCA.pca_tools.default_layout, continuous_cmap=continuous_cmap, discrete_cmap=discrete_cmap
        )
        return fig

    def setup_callbacks(self, app):
        output = Output(self.parent_cid.to_dict(), "figure")

        inputs = dict(
            hue=Input(self.select_hue_cid.to_dict(), "value"),
            pc_x=Input(self.select_pcx_cid.to_dict(), "value"),
            pc_y=Input(self.select_pcy_cid.to_dict(), "value"),
            continuous_cmap=Input(self.select_continuous_cmap_cid.to_dict(), "value"),
            discrete_cmap=Input(self.select_discrete_cmap_cid.to_dict(), "value")
        )
        @app.dash_app.callback(output=output, inputs=inputs)
        def _(hue, pc_x, pc_y, continuous_cmap, discrete_cmap):
            return self.plot(hue, pc_x-1, pc_y-1, continuous_cmap, discrete_cmap)

class PCAProjection(DashPlot):
    def __init__(self, dataset, page_id, loc_class):
        super().__init__(dataset, page_id, loc_class)

        hue_options = self.dataset.get_obs_features(include_genelists=True)

        self.children.update(
            select_continuous_cmap=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-continuous_cmap"),
                options=components.continuous_colormaps, default="viridis",
            ),
            select_discrete_cmap=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-discrete_cmap"),
                options=components.discrete_colormaps, default="scanpy default",
            ),
            select_hue=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-hue"),
                options=hue_options, default=hue_options[0],
                options_callback=lambda: self.dataset.get_obs_features(include_genelists=True),
                update_store_id="update_store-color"
            ),
            input_pcx=InputField(
                cid=CID(self.page_id, self.loc_class, "input-pcx"), type="number",
                default=1, min=1, max=self.dataset.adata.uns["pca"]["variance_ratio"].shape[0],
            ),
            input_pcy=InputField(
                cid=CID(self.page_id, self.loc_class, "input-pcy"), type="number",
                default=2, min=1, max=self.dataset.adata.uns["pca"]["variance_ratio"].shape[0],
            )
        )

        self.actions.update(
            plot_projection=PlotPCA(
                parent_cid=self.cid, dataset=self.dataset,
                select_hue_cid=self.children["select_hue"].cid,
                select_pcx_cid=self.children["input_pcx"].cid,
                select_pcy_cid=self.children["input_pcy"].cid,
                select_discrete_cmap_cid=self.children["select_discrete_cmap"].cid,
                select_continuous_cmap_cid=self.children["select_continuous_cmap"].cid,
            ),
        )

    def create_layout(self) -> list:
        type_params = components.FigureHeaderTab(self.page_id, self.loc_class, tab_label="Type", content=[
            html.Div([
                html.Label("Color"),
                self.children["select_hue"].create_layout(),
                self.children["select_hue"].get_stores()
            ], className="param-row-stacked"),
            # X-axis component
            html.Div([
                html.Label("X Component"),
                self.children["input_pcx"].create_layout(),
                self.children["input_pcx"].get_stores()
            ], className="param-row-stacked"),
            # Y-axis component
            html.Div([
                html.Label("Y Component"),
                self.children["input_pcy"].create_layout(),
                self.children["input_pcy"].get_stores()
            ], className="param-row-stacked"),
        ])

        colormap_tab = components.FigureHeaderTab(self.page_id, self.loc_class, tab_label="Colormap", content=[
            html.Div([
                html.Label("Continuous Color Map"),
                self.children["select_continuous_cmap"].create_layout(),
                self.children["select_continuous_cmap"].get_stores()
            ], className="param-row-stacked"),
            html.Div([
                html.Label("Discrete Color Map"),
                self.children["select_discrete_cmap"].create_layout(),
                self.children["select_discrete_cmap"].get_stores()
            ], className="param-row-stacked"),

        ])

        figure_params = components.FigureHeader(self.page_id, self.loc_class, tabs=[type_params, colormap_tab])

        figure = html.Div([
            html.Div(
                children=figure_params.create_layout(),
                className="fig-header",
            ),
            html.Div([
                dcc.Loading(type="circle", children=[
                    html.Div(
                        dcc.Graph(
                            id=self.cid.to_dict(), className=f"{self.loc_class}-plot"
                        )
                    )
                ])
            ], className=f"{self.loc_class}-body"),
        ], className=f"{self.loc_class}")

        return figure
