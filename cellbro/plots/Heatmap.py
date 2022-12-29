import dash
from dash import Input, Output, State, dcc, html, ALL, ctx
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go

from ..util.Param import Param, ParamsDict
from ..components import components
from ..components.DashFigure import DashFigure
from ..util.DashAction import DashAction
from ..components.DropDown import DropDown

import scout

default_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=5, b=5, l=5, r=5),
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
            allowed_values=components.continuous_colormaps,
        ),
    ]
)

class AddGenesFromList(DashAction):
    def setup_callbacks(self, app):
        output = Output(dict(
            id="input-store", component_id=f"{self.page_id_prefix}-heatmap-selected_genes",
        ), "data")

        inputs = Input(f"{self.page_id_prefix}-heatmap-selected_genelists", "value")

        state = State(f"{self.page_id_prefix}-heatmap-selected_genes", "value")

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(selected_gene_lists, selected_genes):
            if selected_gene_lists is None or selected_genes is None:
                raise PreventUpdate

            if selected_genes is None:
                selected_genes = []

            for gene_list in selected_gene_lists:
                print(self.dataset.uns["gene_lists"])
                selected_genes.extend(self.dataset.get_genes_from_list(gene_list))

            print(selected_gene_lists)
            print(selected_genes)
            raise PreventUpdate
            return dict(value=list(set(selected_genes)))

class EditGeneList(DashAction):
    def setup_callbacks(self, app):
        @app.dash_app.callback(
            output=Output(dict(
                id="input-store", component_id=f"{self.page_id_prefix}-heatmap-selected_genelists",
            ), "data"),

            inputs=Input("genelist-store", "data")
        )
        def _(genelist_store):
            return dict(options=self.dataset.get_gene_lists())

class PlotHeatmap(DashAction):
    def plot(self, params):
        selected_genes = params.pop("selected_genes")

        if selected_genes and len(selected_genes) > 0:
            selected_genes = selected_genes
        else:
            selected_genes = self.dataset.adata.var_names[:50].tolist()

        fig = scout.ply.heatmap(
            adata=self.dataset.adata, var_names=selected_genes, categoricals=params["categoricals"],
            layer=params["layer"], cluster_cells_by=params["cluster_cells_by"], layout=dict(),
            cmap=params["colormap"],
        )

        style = {
            # "width": f"{int(z.shape[0]/2)+100}px",
            "height": fig.layout["height"],
        }
        return fig, style

    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-heatmap-plot", component_property="figure"),
            Output(component_id=f"{self.page_id_prefix}-heatmap-plot", component_property="style"),
            Output(f"{self.page_id_prefix}-heatmap-submit", "n_clicks"),
        ]

        # Inputs to Projection
        inputs = dict(
            submit=Input(f"{self.page_id_prefix}-heatmap-submit", "n_clicks"),
            store=Input(dict(id="selected-store", component_id=ALL), "data"),
        )

        state = dict(
            selected_genes=State(f"{self.page_id_prefix}-heatmap-selected_genes", "value"),
            cluster_cells_by=State(f"{self.page_id_prefix}-heatmap-cluster_cells_by", "value"),
            categoricals=State(f"{self.page_id_prefix}-heatmap-selected_categoricals", "value"),
        )

        for key in heatmap_params.keys():
            state[key] = State(component_id=f"{self.page_id_prefix}-heatmap-{key}", component_property="value")
            

        @app.dash_app.callback(output=output, inputs=inputs, state=state)
        def _(submit, store, **kwargs):
            if submit is None or ctx.triggered_id == f"{self.page_id_prefix}-heatmap-submit":
                if kwargs["selected_genes"] is None or len(kwargs["selected_genes"]) == 0:
                    return dash.no_update, dash.no_update, 0

                fig, style = self.plot(kwargs)
                return fig, style, 0

            raise PreventUpdate

class Heatmap(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)

        genes = sorted(self.dataset.adata.var_names.tolist())
        gene_lists = sorted(self.dataset.get_gene_lists())
        categoricals = self.dataset.get_categoric()

        self.select_genes = DropDown(
            self.page_id_prefix, id=f"{self.page_id_prefix}-heatmap-selected_genes",
            options=genes, default=None, clearable=True,
            placeholder="Select Genes", multi=True,
            style={"width": "100%"}
        )
        self.select_genelist = DropDown(
            self.page_id_prefix, id=f"{self.page_id_prefix}-heatmap-selected_genelists",
            clearable=True, default=None, options=gene_lists,
            placeholder="Select Gene Lists", multi=True,
            style={"width": "100%"}
        )
        self.select_clusterby = DropDown(
            self.page_id_prefix, id=f"{self.page_id_prefix}-heatmap-cluster_cells_by",
            options=categoricals + ["barcode"], default=None, clearable=True,
            placeholder="Select Categorical (Optional)",
        )
        self.select_categoricals = DropDown(
            self.page_id_prefix, id=f"{self.page_id_prefix}-heatmap-selected_categoricals",
            options=categoricals, default=next(iter(categoricals), None),
            placeholder="Select Feature(s)", multi=True, clearable=True,
            style={"width": "100%"}
        )
        self.select_layer = DropDown(
            self.page_id_prefix, id=f"{self.page_id_prefix}-heatmap-layer",
            options=heatmap_params["layer"].allowed_values, default=heatmap_params["layer"].default,
        )
        self.select_colormap = DropDown(
            self.page_id_prefix, id=f"{self.page_id_prefix}-heatmap-colormap",
            options=heatmap_params["colormap"].allowed_values, default=heatmap_params["colormap"].default,
        )

        self.actions.update(self.select_genes.actions)
        self.actions.update(self.select_genelist.actions)
        self.actions.update(self.select_clusterby.actions)
        self.actions.update(self.select_categoricals.actions)
        self.actions.update(self.select_layer.actions)
        self.actions.update(self.select_colormap.actions)
        self.actions.update(
            plot_heatmap=PlotHeatmap(dataset, self.page_id_prefix, self.loc_class),
            add_genes_from_list=AddGenesFromList(dataset, self.page_id_prefix, self.loc_class),
            edit_gene_list=EditGeneList(dataset, self.page_id_prefix, self.loc_class),
        )

    def get_sidebar_params(self) -> list:
        divs = [
            html.Div([
                html.Label("Show Genes",className="param-label"),
                html.Div([
                    self.select_genes.create_layout(),
                    self.select_genes.get_stores(),
                    self.select_genelist.create_layout(),
                    self.select_genelist.get_stores(),
                ], style={"display": "flex", "gap": "10px"})
            ], className="param-row-stacked"),
            html.Div(
                children=[
                    html.Label("Cluster Group By", className="param-label"),
                    html.Div([
                        self.select_clusterby.create_layout(),
                        self.select_clusterby.get_stores(),
                    ], className="param-select"),
                ],
                className="param-row-stacked",
            ),
            html.Div(
                children=[
                    html.Label("Show Categorical Feature(s)", className="param-label"),
                    html.Div([
                        self.select_categoricals.create_layout(),
                        self.select_categoricals.get_stores(),
                    ], className="param-select"),
                ],
                className="param-row-stacked",
            ),
            html.Div(
                children=[
                    html.Label("Select Layer", className="param-label"),
                    html.Div([
                        self.select_layer.create_layout(),
                        self.select_layer.get_stores(),
                    ], className="param-select"),
                ],
                className="param-row-stacked",
            ),
            html.Div(
                children=[
                    html.Label("Select Colormap", className="param-label"),
                    html.Div([
                        self.select_colormap.create_layout(),
                        self.select_colormap.get_stores(),
                    ], className="param-select"),
                ],
                className="param-row-stacked",
            ),
        ]
        # divs.extend(components.params_layout(heatmap_params, f"{self.page_id_prefix}-heatmap"))


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
                    className=f"{self.loc_class}-body", id=f"heatmap-figure"
                )
            ]
        )

        return figure_layout
