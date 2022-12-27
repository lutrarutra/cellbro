import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

from cellbro.util.Param import *
import cellbro.util.Components as Components
from cellbro.util.DashAction import DashAction
from cellbro.plots.DashFigure import DashFigure

import scout

heatmap_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=10, b=10, l=10, r=10),
)

heatmap_params = ParamsDict(
    [
        Param(
            key="layer",
            name="Layer",
            default="logcentered",
            type=list,
            description="",
            allowed_values={
                "log1p": "log1p",
                "counts": "Counts",
                "ncounts": "Normalized Counts",
                "centered": "Centered Counts",
                "logcentered": "Log Centered",
            },
        ),
        Param(
            key="colormap",
            name="Colormap",
            default="seismic",
            type=list,
            description="",
            allowed_values=Components.continuous_colormaps,
        ),
    ]
)

class AddGenesFromList(DashAction):
    def setup_callbacks(self, app):
        output = Output(f"{self.page_id_prefix}-heatmap-selected_genes", "value")
        inputs = [Input(f"{self.page_id_prefix}-heatmap-selected_genelists", "value")]
        state = [State(f"{self.page_id_prefix}-heatmap-selected_genes", "value")]
        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(gene_lists, selected_genes):
            if gene_lists is None:
                raise PreventUpdate

            if selected_genes is None:
                selected_genes = []

            for gene_list in gene_lists:
                selected_genes.extend(self.dataset.get_genes_from_list(gene_list))

            return list(set(selected_genes))

class PlotHeatmap(DashAction):
    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-heatmap-plot", component_property="figure"),
            Output(component_id=f"{self.page_id_prefix}-heatmap-plot", component_property="style"),
        ]

        # Inputs to Projection
        inputs = {
            "submit": Input(
                component_id=f"{self.page_id_prefix}-heatmap-submit", component_property="n_clicks"
            ),
        }

        state = dict(
            selected_genes=State(f"{self.page_id_prefix}-heatmap-selected_genes", "value"),
            cluster_cells_by=State(f"{self.page_id_prefix}-heatmap-cluster_cells_by", "value"),
            categoricals=State(f"{self.page_id_prefix}-heatmap-selected_categoricals", "value"),
        )

        for key in heatmap_params.keys():
            state[key] = State(component_id=f"{self.page_id_prefix}-heatmap-{key}", component_property="value")
            

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, **kwargs):
            return Heatmap.plot(self.dataset, kwargs)

class EditGeneList(DashAction):
    def setup_callbacks(self, app):
        @app.dash_app.callback(
            output=Output(f"{self.page_id_prefix}-heatmap-selected_genelists", "options"),
            inputs=[Input("genelist-store", "data")]
        )
        def _(genelist_store):
            return self.dataset.get_gene_lists()


class Heatmap(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.actions.update(
            plot_heatmap=PlotHeatmap(dataset, self.page_id_prefix),
            add_genes_from_list=AddGenesFromList(dataset, self.page_id_prefix),
            edit_gene_list=EditGeneList(dataset, self.page_id_prefix),
        )

    @staticmethod
    def plot(dataset, params):
        selected_genes = params.pop("selected_genes")

        if selected_genes and len(selected_genes) > 0:
            selected_genes = selected_genes
        else:
            selected_genes = dataset.adata.var_names[:50].tolist()

        fig = scout.ply.heatmap(
            adata=dataset.adata, var_names=selected_genes, categoricals=params["categoricals"],
            layer=params["layer"], cluster_cells_by=params["cluster_cells_by"], layout=dict()
        )

        style = {
            # "width": f"{int(z.shape[0]/2)+100}px",
            "height": fig.layout["height"],
        }
        return fig, style

    def get_sidebar_params(self) -> list:
        genes = sorted(self.dataset.adata.var_names.tolist())
        gene_lists = sorted(self.dataset.get_gene_lists())
        categoricals = self.dataset.get_categoric()

        divs = [
            html.Div([
                html.Label(
                    "Show Genes",
                    className="param-label",
                    ),
                html.Div([
                    dcc.Dropdown(
                        options=genes, value=None, id=f"{self.page_id_prefix}-heatmap-selected_genes", clearable=True,
                        placeholder="Select Genes", multi=True,
                        style={"width": "100%"}
                        ),
                    dcc.Dropdown(
                        options=gene_lists, value=None, id=f"{self.page_id_prefix}-heatmap-selected_genelists", clearable=True,
                        placeholder="Select Gene Lists", multi=True,
                        style={"width": "100%"}
                    ),
                ], style={"display": "flex", "gap": "10px"})
            ], className="param-row-stacked"),
            html.Div(
                children=[
                    html.Label(
                        "Cluster Group By",
                        className="param-label",
                    ),
                    html.Div(
                        [
                            dcc.Dropdown(
                                options=categoricals + ["barcode"], value=None,
                                id=f"{self.page_id_prefix}-heatmap-cluster_cells_by", clearable=True,
                                placeholder="Select Categorical (Optional)",
                            )
                        ],
                        className="param-select",
                    ),
                ],
                className="param-row-stacked",
            ),
            html.Div(
                children=[
                    html.Label(
                        "Show Categorical Feature(s)",
                        className="param-label",
                    ),
                    html.Div(
                        [
                            dcc.Dropdown(
                                options=categoricals, value=next(iter(categoricals), None),
                                id=f"{self.page_id_prefix}-heatmap-selected_categoricals",
                                placeholder="Select Feature(s)", multi=True, clearable=True,
                                style={"width": "100%"}
                            ),
                        ],
                        className="param-select",
                    ),
                ],
                className="param-row-stacked",
            ),
        ]
        divs.extend(Components.params_layout(heatmap_params, f"{self.page_id_prefix}-heatmap"))

        return divs

    def create_layout(self):
        figure_layout = html.Div(
            children=[
                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    [
                                        dcc.Graph(
                                            id=f"{self.page_id_prefix}-heatmap-plot", className=f"{self.loc_class}-plot"
                                        )
                                    ],
                                )
                            ],
                        )
                    ],
                    className=f"{self.loc_class}-body",
                )
            ]
        )

        return figure_layout
