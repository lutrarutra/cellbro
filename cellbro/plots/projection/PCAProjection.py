import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

from ...components.DashPlot import DashPlot
from ...util.DashAction import DashAction
from ...components import components
from .. import PCA

import scout

class PlotPCA(DashAction):
    def plot(self, color, pc_x, pc_y, continuous_cmap, discrete_cmap):
        fig = scout.ply.projection(
            self.dataset.adata, obsm_layer="X_pca", hue=color, components=[pc_x, pc_y],
            layout=PCA.pca_tools.default_layout, continuous_cmap=continuous_cmap, discrete_cmap=discrete_cmap
        )
        return fig

    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id}-main-plot", component_property="figure"),
        ]

        inputs = {
            "color": Input(component_id=f"{self.page_id}-projection-color", component_property="value"),
            "pc_x": Input(component_id=f"{self.page_id}-projection-x-component", component_property="value"),
            "pc_y": Input(component_id=f"{self.page_id}-projection-y-component", component_property="value"),
            "continuous_cmap": Input(
                component_id=f"{self.page_id}-projection-continuous_cmap", component_property="value"
            ),
            "discrete_cmap": Input(
                component_id=f"{self.page_id}-projection-discrete_cmap", component_property="value"
            )
        }

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(color, pc_x, pc_y, continuous_cmap, discrete_cmap):
            return [self.plot(color, pc_x-1, pc_y-1, continuous_cmap, discrete_cmap)]

class PCAProjection(DashPlot):
    def __init__(self, dataset, page_id, loc_class):
        super().__init__(dataset, page_id, loc_class)
        self.actions.update(
            plot_projection=PlotPCA(self.dataset, self.page_id, loc_class),
        )

    def create_layout(self) -> list:
        type_params = components.FigureHeaderTab(self.page_id, tab_label="Type", children=[
            html.Div([
                html.Label("Color"),
                dcc.Dropdown(
                    self.dataset.adata.obs_keys() + list(self.dataset.adata.var_names),
                    value=self.dataset.adata.obs_keys()[0],
                    id=f"{self.page_id}-projection-color", clearable=False,
                ),
            ], className="param-row-stacked"),
            # X-axis component
            html.Div([
                html.Label("X Component"),
                dbc.Input(
                    id=f"{self.page_id}-projection-x-component", type="number", value=1, min=1, step=1,
                    max=self.dataset.adata.uns["pca"]["variance_ratio"].shape[0] + 1,
                    className="param-input",
                ),
            ], className="param-row-stacked"),
            # X-axis component
            html.Div(
                children=[
                    html.Label("Y Component"),
                    dbc.Input(
                        id=f"{self.page_id}-projection-y-component", type="number", value=2, min=1, step=1,
                        max=self.dataset.adata.uns["pca"]["variance_ratio"].shape[0] + 1,
                        className="param-input",
                    ),
                ],
                className="param-row-stacked",
            ),
        ])

        colormap_tab = components.FigureHeaderTab(self.page_id, tab_label="Colormap", children=[
            html.Div([
                html.Label("Continuous Color Map"),
                components.create_colormap_selector(
                    id=f"{self.page_id}-projection-continuous_cmap",
                    options=components.continuous_colormaps,
                    default="viridis",
                )
            ], className="param-row-stacked"),
            html.Div([
                html.Label("Discrete Color Map"),
                components.create_colormap_selector(
                    id=f"{self.page_id}-projection-discrete_cmap",
                    options=components.discrete_colormaps,
                    default="scanpy default",
                )
            ], className="param-row-stacked"),

        ])

        figure_params = components.FigureHeader(self.page_id, tabs=[type_params, colormap_tab])

        figure = html.Div(
            children=[
                html.Div(
                    children=figure_params.create_layout(),
                    className="fig-header",
                ),
                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id=f"{self.page_id}-{self.loc_class}-plot",
                                        className=f"{self.loc_class}-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    className=f"{self.loc_class}-body",
                ),
            ],
            className=f"{self.loc_class}",
        )
        return figure
