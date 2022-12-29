import dash_bootstrap_components as dbc
import plotly.express as px
import scanpy as sc
from dash import Input, Output, State, dcc, html, ctx
from dash.exceptions import PreventUpdate

from ...components.DashFigure import DashFigure
from ...components import components
from ...util.DashAction import DashAction
from ...components.DropDown import DropDown

from .UMAP import UMAP, SCVI_UMAP
from .TSNE import TSNE
from .Trimap import Trimap
from . import prj_tools

import scout

# class UpdateColorOptions(DashAction):
#     def setup_callbacks(self, app):
#         output = [
#             Output(f"{self.page_id_prefix}-projection-color", "options"),
#             Output(f"{self.page_id_prefix}-projection-color", "value"),
#         ]

#         inputs = dict(
#             store=Input("genelist-store", "data")
#         )
#         @app.dash_app.callback(output, inputs)
#         def _(store):
#             if store is None:
#                 raise PreventUpdate

#             color_options = self.dataset.get_obs_features(include_gene_lists=True)
#             return color_options, next(iter(color_options), None)


class ApplyProjection(DashAction):
    def apply(self, projection_type, params):
        if projection_type == "UMAP":
            projection = UMAP(self.dataset, UMAP.parse_params(params))
        elif projection_type == "t-SNE":
            projection = TSNE(self.dataset, TSNE.parse_params(params))
        elif projection_type == "Trimap":
            projection = Trimap(self.dataset, Trimap.parse_params(params))

        return projection.apply()

    def setup_callbacks(self, app):
        output = Output(dict(
            id=f"input-store",
            component_id=f"{self.page_id_prefix}-projection-plot-type"
        ), "data")

        inputs = dict(projection_submit=Input(f"{self.page_id_prefix}-projection-submit", "n_clicks"))

        state = dict()
        state["projection_type"] = State(f"{self.page_id_prefix}-projection-type-select", "value")
        for key in UMAP._params.keys():
            state[f"umap_{key}"] = State(f"{self.page_id_prefix}-projection-umap-{key}", "value")

        for key in SCVI_UMAP._params.keys():
            state[f"scvi_umap_{key}"] = State(f"{self.page_id_prefix}-projection-scvi_umap-{key}", "value")

        for key in TSNE._params.keys():
            state[f"tsne_{key}"] = State(f"{self.page_id_prefix}-projection-tsne-{key}", "value")

        for key in Trimap._params.keys():
            state[f"trimap_{key}"] = State(f"{self.page_id_prefix}-projection-trimap-{key}", "value")

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(projection_submit, projection_type, **kwargs):
            if projection_submit is None:
                raise PreventUpdate

            value = self.apply(projection_type, params=kwargs)
            options = list(self.dataset.adata.obsm.keys())
            return dict(value=value, options=options)

class PlotProjection(DashAction):
    def plot(self, color, obsm_layer, continuous_cmap, discrete_cmap):
        if "Gene List:" in color:
            gene_list = color.split("Gene List: ")[1]
            color = self.dataset.adata.uns["gene_lists"][gene_list]

        fig = scout.ply.projection(
            self.dataset.adata, obsm_layer=obsm_layer, hue=color,
            layout=prj_tools.default_layout, continuous_cmap=continuous_cmap, discrete_cmap=discrete_cmap
        )
        return fig

    def setup_callbacks(self, app):
        outputs = Output(component_id=f"{self.page_id_prefix}-projection-plot", component_property="figure")

        inputs = dict(
            color=Input(f"{self.page_id_prefix}-projection-color", "value"),
            obsm_layer=Input(f"{self.page_id_prefix}-projection-plot-type", "value"),
            continuous_cmap=Input(f"{self.page_id_prefix}-projection-continuous_cmap", "value"),
            discrete_cmap=Input(f"{self.page_id_prefix}-projection-discrete_cmap", "value"),
        )
        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(color, obsm_layer, continuous_cmap, discrete_cmap):
            if obsm_layer not in self.dataset.adata.obsm.keys():
                raise PreventUpdate
            fig = self.plot(
                color=color, obsm_layer=obsm_layer, continuous_cmap=continuous_cmap,
                discrete_cmap=discrete_cmap
            )
            return fig

class SelectProjectionType(DashAction):
    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-projection-umap", component_property="style"),
            Output(component_id=f"{self.page_id_prefix}-projection-tsne", component_property="style"),
            Output(component_id=f"{self.page_id_prefix}-projection-trimap", component_property="style"),
        ]
        inputs = [
            Input(component_id=f"{self.page_id_prefix}-projection-type-select", component_property="value"),
        ]

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(projection_type):
            if projection_type == "UMAP":
                return {"display": "block"}, {"display": "none"}, {"display": "none"}

            if projection_type == "t-SNE":
                return {"display": "none"}, {"display": "block"}, {"display": "none"}

            if projection_type == "Trimap":
                return {"display": "none"}, {"display": "none"}, {"display": "block"}


class Projection(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.actions.update(
            plot_projection=PlotProjection(self.dataset, self.page_id_prefix, "main"),
            apply_projection=ApplyProjection(self.dataset, self.page_id_prefix, "main"),
            # update_color_options=UpdateColorOptions(self.dataset, self.page_id_prefix, self.loc_class),
            select_projection_type=SelectProjectionType(self.dataset, self.page_id_prefix, self.loc_class),
        )

        available_projections = list(self.dataset.adata.obsm.keys())
        color_options = self.dataset.get_obs_features(include_gene_lists=True)

        self.select_projection_type = DropDown(
            self.page_id_prefix, id=f"{self.page_id_prefix}-projection-plot-type",
            options=available_projections, default=available_projections[0]
        )
        self.select_color = DropDown(
            self.page_id_prefix, id=f"{self.page_id_prefix}-projection-color",
            options=color_options, default=color_options[0],
        )

        self.select_continuous_cmap = DropDown(
            self.page_id_prefix, id=f"{self.page_id_prefix}-projection-continuous_cmap",
            options=components.continuous_colormaps, default="viridis",
        )

        self.select_discrete_cmap = DropDown(
            self.page_id_prefix, id=f"{self.page_id_prefix}-projection-discrete_cmap",
            options=components.discrete_colormaps, default="scanpy default",
        )

        self.actions.update(self.select_projection_type.actions)
        self.actions.update(self.select_color.actions)
        self.actions.update(self.select_continuous_cmap.actions)
        self.actions.update(self.select_discrete_cmap.actions)
        # print(self.actions.keys())


    def create_layout(self) -> list:
        projection_type_tab = components.FigureHeaderTab(self.page_id_prefix, tab_label="Type", children=[
            # Projection type celect
            html.Div([
                html.Label("Projection Type"),
                self.select_projection_type.create_layout(),
                self.select_projection_type.get_stores(),
                # dcc.Dropdown(
                #     available_projections, value=available_projections[0],
                #     id=f"{self.page_id_prefix}-projection-plot-type", clearable=False,
                #     persistence=True, persistence_type="local",
                # ),
            ], className="param-row-stacked"),
            # Projection Hue celect
            html.Div([
                html.Label("Color"),
                self.select_color.create_layout(),
                self.select_color.get_stores(),
                # dcc.Dropdown(
                #     options=color_options,
                #     value=color_options[0],
                #     id=f"{self.page_id_prefix}-projection-color",
                #     clearable=False,
                # )
            ], className="param-row-stacked")
        ])

        colormap_tab = components.FigureHeaderTab(self.page_id_prefix, tab_label="Colormap", children=[
            html.Div([
                html.Label("Continuous Color Map"),
                self.select_continuous_cmap.create_layout(),
                self.select_continuous_cmap.get_stores(),
                # components.create_colormap_selector(
                #     id=f"{self.page_id_prefix}-projection-continuous_cmap",
                #     options=components.continuous_colormaps,
                #     default="viridis",
                # )
            ], className="param-row-stacked"),
            html.Div([
                html.Label("Discrete Color Map"),
                self.select_discrete_cmap.create_layout(),
                self.select_discrete_cmap.get_stores(),
                # components.create_colormap_selector(
                #     id=f"{self.page_id_prefix}-projection-discrete_cmap",
                #     options=components.discrete_colormaps,
                #     default="scanpy default",
                # )
            ], className="param-row-stacked")
        ])

        fig_params = components.FigureHeader(self.page_id_prefix, tabs=[projection_type_tab, colormap_tab])

        figure_layout = html.Div(
            children=[
                html.Div(children=fig_params.create_layout(), className="fig-header"),

                html.Div([
                    dcc.Loading(type="circle", children=[
                        html.Div(
                            dcc.Graph(
                                id=f"{self.page_id_prefix}-projection-plot", className=f"{self.loc_class}-plot"
                            )
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
                            id=f"{self.page_id_prefix}-projection-type-select",
                            options=["UMAP", "t-SNE", "Trimap"],
                            value="UMAP", clearable=False,
                        )
                    ],
                    className="param-row-stacked",
                )
            ]),
            html.Div(
                children=UMAP.get_layout(self.page_id_prefix),
                style={"display": "none"},
                id=f"{self.page_id_prefix}-projection-{UMAP.get_key()}",
                className="param-class"
            ),
            html.Div(
                children=TSNE.get_layout(self.page_id_prefix),
                style={"display": "none"},
                id=f"{self.page_id_prefix}-projection-{TSNE.get_key()}",
                className="param-class"
            ),
            html.Div(
                children=Trimap.get_layout(self.page_id_prefix),
                style={"display": "none"},
                id=f"{self.page_id_prefix}-projection-{Trimap.get_key()}",
                className="param-class"
            ),
            html.Div(
                children=SCVI_UMAP.get_layout(self.page_id_prefix),
                style={"display": "none"},
                id=f"{self.page_id_prefix}-projection-{SCVI_UMAP.get_key()}",
                className="param-class"
            ),
        ]

