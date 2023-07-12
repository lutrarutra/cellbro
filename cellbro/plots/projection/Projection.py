import dash_bootstrap_components as dbc
import plotly.express as px
import scanpy as sc
from dash import dcc, html, ctx, Input, Output, State, ALL
from dash.exceptions import PreventUpdate

from ...components.DashPlot import DashPlot
from ...components import components
from ...util.DashAction import DashAction
from ...components.DropDown import DropDown
from ...components.DashComponent import DashComponent
from ...components.CID import CID, LocClass

from .UMAP import UMAP, SCVI_UMAP
from .TSNE import TSNE
from .Trimap import Trimap
from . import prj_tools

import scout

class ApplyProjection(DashAction):
    def __init__(
        self, parent_cid: CID, dataset
    ):
        super().__init__(parent_cid, dataset)

    def apply(self, projection_type, params):
        if projection_type == "UMAP":
            projection = UMAP(self.dataset, UMAP.parse_params(params))
        elif projection_type == "t-SNE":
            projection = TSNE(self.dataset, TSNE.parse_params(params))
        elif projection_type == "Trimap":
            projection = Trimap(self.dataset, Trimap.parse_params(params))

        return projection.apply()

    def setup_callbacks(self, app):
        output = [
            Output("update_store-projection_type", "data"),
            Output(dict(page_id=self.page_id, loc_class=self.loc_class.name, type="select-projection_type"), "value"),
        ]

        inputs = dict(
            projection_submit=Input(f"{self.page_id}-{self.loc_class}-sidebar-apply_btn", "n_clicks"),
            projection_type=State(f"{self.page_id}-projection-type-select", "value"),
        )

        for key in UMAP._params.keys():
            inputs[f"umap_{key}"] = State(f"{self.page_id}-projection-umap-{key}", "value")

        for key in SCVI_UMAP._params.keys():
            inputs[f"scvi_umap_{key}"] = State(f"{self.page_id}-projection-scvi_umap-{key}", "value")

        for key in TSNE._params.keys():
            inputs[f"tsne_{key}"] = State(f"{self.page_id}-projection-tsne-{key}", "value")

        for key in Trimap._params.keys():
            inputs[f"trimap_{key}"] = State(f"{self.page_id}-projection-trimap-{key}", "value")

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(projection_submit, projection_type, **kwargs):
            if projection_submit is None:
                raise PreventUpdate

            value = self.apply(projection_type, params=kwargs)
            options = list(self.dataset.adata.obsm.keys())
            return dict(options=options), value

class SelectProjectionType(DashAction):
    def setup_callbacks(self, app):
        output = [
            Output(f"{self.page_id}-projection-umap", "style"),
            Output(f"{self.page_id}-projection-tsne", "style"),
            Output(f"{self.page_id}-projection-trimap", "style"),
        ]
        inputs = [
            Input(f"{self.page_id}-projection-type-select", "value"),
        ]

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(projection_type):
            if projection_type == "UMAP":
                return {"display": "block"}, {"display": "none"}, {"display": "none"}

            if projection_type == "t-SNE":
                return {"display": "none"}, {"display": "block"}, {"display": "none"}

            if projection_type == "Trimap":
                return {"display": "none"}, {"display": "none"}, {"display": "block"}

class PlotProjection(DashAction):
    def __init__(
        self, parent_cid: CID, dataset,
        select_color_cid: CID,
        obsm_layer_cid: CID,
        continuous_cmap_cid: CID,
        discrete_cmap_cid: CID,
    ):
        super().__init__(parent_cid, dataset)
        self.select_color_cid = select_color_cid
        self.obsm_layer_cid = obsm_layer_cid
        self.continuous_cmap_cid = continuous_cmap_cid
        self.discrete_cmap_cid = discrete_cmap_cid

    def plot(self, color, obsm_layer, continuous_cmap, discrete_cmap):
        if "Gene List:" in color:
            genelist = color.split("Gene List: ")[1]
            color = self.dataset.adata.uns["genelists"][genelist]

        fig = scout.ply.projection(
            self.dataset.adata, obsm_layer=obsm_layer, color=color,
            layout=prj_tools.default_layout, continuous_cmap=continuous_cmap, discrete_cmap=discrete_cmap
        )
        return fig

    def setup_callbacks(self, app):
        @app.dash_app.callback(
            output=Output(self.parent_cid.to_dict(), "figure"),
            inputs=[
                Input(self.select_color_cid.to_dict(), "value"),
                Input(self.obsm_layer_cid.to_dict(), "value"),
                Input(self.continuous_cmap_cid.to_dict(), "value"),
                Input(self.discrete_cmap_cid.to_dict(), "value")
            ]
        )
        def _(color, obsm_layer, continuous_cmap, discrete_cmap):
            if obsm_layer not in self.dataset.adata.obsm.keys():
                raise PreventUpdate

            fig = self.plot(
                color=color, obsm_layer=obsm_layer, continuous_cmap=continuous_cmap,
                discrete_cmap=discrete_cmap
            )
            return fig

class Projection(DashPlot):
    def __init__(self, dataset, page_id, loc_class):
        super().__init__(dataset, page_id, loc_class)

        available_projections = list(self.dataset.adata.obsm.keys())
        color_options = self.dataset.get_obs_features(include_genelists=True)

        self.children.update(
            select_projection_type=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-projection_type"),
                options=available_projections, default=available_projections[0],
                options_callback= lambda: list(self.dataset.adata.obsm.keys())
            ),
            select_color=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-color"),
                options=color_options, default=color_options[0],
                options_callback=lambda: self.dataset.get_obs_features(include_genelists=True)
            ),
            select_continuous_cmap=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-continuous_cmap"),
                options=components.continuous_colormaps, default="viridis",
            ),
            select_discrete_cmap=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-discrete_cmap"),
                options=components.discrete_colormaps, default="scanpy default",
            )
        )
        self.children.update(
            type_header_tab=components.FigureHeaderTab(self.page_id, self.loc_class, tab_label="Type", content=[
                # Projection type celect
                html.Div([
                    html.Label("Projection Type"),
                    self.children["select_projection_type"].create_layout(),
                    self.children["select_projection_type"].get_stores(),
                ], className="param-row-stacked"),
                # Projection Hue celect
                html.Div([
                    html.Label("Color"),
                    self.children["select_color"].create_layout(),
                    self.children["select_color"].get_stores(),
                ], className="param-row-stacked")
            ]),
            color_header_tab=components.FigureHeaderTab(self.page_id, self.loc_class, tab_label="Colormap", content=[
                html.Div([
                    html.Label("Continuous Color Map"),
                    self.children["select_continuous_cmap"].create_layout(),
                    self.children["select_continuous_cmap"].get_stores(),
                ], className="param-row-stacked"),
                html.Div([
                    html.Label("Discrete Color Map"),
                    self.children["select_discrete_cmap"].create_layout(),
                    self.children["select_discrete_cmap"].get_stores(),
                ], className="param-row-stacked")
            ])
        )

        self.actions.update(
            plot_projection=PlotProjection(
                parent_cid=self.cid,
                dataset=self.dataset,
                select_color_cid=self.children["select_color"].cid,
                obsm_layer_cid=self.children["select_projection_type"].cid,
                continuous_cmap_cid=self.children["select_continuous_cmap"].cid,
                discrete_cmap_cid=self.children["select_discrete_cmap"].cid,
            ),
            apply_projection=ApplyProjection(
                parent_cid=self.cid, dataset=self.dataset,
            ),
            select_projection_type=SelectProjectionType(
                parent_cid=self.cid, dataset=self.dataset,
            )
        )

    def create_layout(self) -> list:
        fig_params = components.FigureHeader(self.page_id, self.loc_class, tabs=[
            self.children["type_header_tab"],
            self.children["color_header_tab"]
        ])

        figure_layout = html.Div(
            children=[
                html.Div(children=fig_params.create_layout(), className="fig-header"),

                html.Div([
                    dcc.Loading(type="circle", children=[
                        html.Div(
                            dcc.Graph(id=self.cid.to_dict(), className=f"{self.loc_class}-plot")
                        )],
                    )
                ], className=f"{self.loc_class}-body"),
            ],
            className=f"{self.loc_class}",
        )

        return figure_layout

    def get_sidebar_params(self):
        return [
            html.Div([
                html.Div(
                    children=[
                        html.Label(
                            "Projection Type",
                            className="param-label",
                        ),
                        dcc.Dropdown(
                            id=f"{self.page_id}-projection-type-select",
                            options=["UMAP", "t-SNE", "Trimap"],
                            value="UMAP", clearable=False,
                        )
                    ],
                    className="param-row-stacked",
                )
            ]),
            html.Div(
                children=UMAP.get_layout(self.page_id),
                style={"display": "none"},
                id=f"{self.page_id}-projection-{UMAP.get_key()}",
                className="param-class"
            ),
            html.Div(
                children=TSNE.get_layout(self.page_id),
                style={"display": "none"},
                id=f"{self.page_id}-projection-{TSNE.get_key()}",
                className="param-class"
            ),
            html.Div(
                children=Trimap.get_layout(self.page_id),
                style={"display": "none"},
                id=f"{self.page_id}-projection-{Trimap.get_key()}",
                className="param-class"
            ),
            html.Div(
                children=SCVI_UMAP.get_layout(self.page_id),
                style={"display": "none"},
                id=f"{self.page_id}-projection-{SCVI_UMAP.get_key()}",
                className="param-class"
            ),
        ]

