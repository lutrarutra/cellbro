from functools import namedtuple

import dash
from dash import Input, Output, State, dcc, html, ALL, ctx
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go

from ..util.Param import Param, ParamsDict
from ..components import components
from ..components.DashPlot import DashPlot
from ..util.DashAction import DashAction
from ..components.DropDown import DropDown
from ..components.CID import CID

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

# class AddGenesFromList(DashAction):
#     def setup_callbacks(self, app):
#         output = Output(dict(
#             id="input-store", component_id=f"{self.page_id}-{self.loc_class}-param-genes",
#         ), "data")

#         inputs = Input(f"{self.page_id}-{self.loc_class}-param-genelists", "value")

#         state = State(f"{self.page_id}-{self.loc_class}-param-genes", "value")

#         @app.dash_app.callback(output=output, inputs=inputs, state=state)
#         def _(selected_genelists, selected_genes):
#             if selected_genelists is None or selected_genes is None:
#                 raise PreventUpdate

#             if selected_genes is None:
#                 selected_genes = []

#             for genelist in selected_genelists:
#                 selected_genes.extend(self.dataset.get_genes_from_list(genelist))

#             raise PreventUpdate
#             return dict(value=list(set(selected_genes)))

class EditGeneList(DashAction):
    def setup_callbacks(self, app):
        @app.dash_app.callback(
            output=Output(dict(
                id="input-store", component_id=f"{self.page_id}-heatmap-selected_genelists",
            ), "data"),

            inputs=Input("genelist-store", "data")
        )
        def _(genelist_store):
            return dict(options=self.dataset.get_genelists())

class PlotHeatmap(DashAction):
    RType = namedtuple("RType", ["figure", "style", "button"])

    def __init__(
        self, parent_cid: CID, dataset,
        select_genes_cid,
        select_clusterby_cid,
        select_categoricals_cid,
        select_layer_cid,
        select_colormap_cid,
    ):
        super().__init__(parent_cid, dataset)
        self.select_genes_cid = select_genes_cid
        self.select_clusterby_cid = select_clusterby_cid
        self.select_categoricals_cid = select_categoricals_cid
        self.select_layer_cid = select_layer_cid
        self.select_colormap_cid = select_colormap_cid

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
            Output(f"{self.page_id}-{self.loc_class}-sidebar-apply_btn", "n_clicks"),
        ]

        # Inputs to Projection
        inputs = dict(
            submit=Input(f"{self.page_id}-{self.loc_class}-sidebar-apply_btn", "n_clicks"),
            selected_genes=Input(self.select_genes_cid.to_dict(), "value"),
            cluster_cells_by=Input(self.select_clusterby_cid.to_dict(), "value"),
            categoricals=Input(self.select_categoricals_cid.to_dict(), "value"),
            layer=Input(self.select_layer_cid.to_dict(), "value"),
            colormap=Input(self.select_colormap_cid.to_dict(), "value"),
        )
            

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(submit, selected_genes, cluster_cells_by, categoricals, layer, colormap):
            if submit is None or ctx.triggered_id == f"{self.page_id}-{self.loc_class}-sidebar-apply_btn":
                if (
                    selected_genes is None or len(selected_genes) == 0
                    ):
                    return self.RType(
                        figure=dash.no_update, style=dash.no_update, button=None
                    )

                fig, style = self.plot(selected_genes, cluster_cells_by, categoricals, layer, colormap)

                return self.RType(
                    figure=fig, style=style, button=0
                )

            raise PreventUpdate

class Heatmap(DashPlot):
    def __init__(self, dataset, page_id, loc_class):
        super().__init__(dataset, page_id, loc_class)

        genes = sorted(self.dataset.adata.var_names.tolist())
        genelists = sorted(self.dataset.get_genelists())
        categoricals = self.dataset.get_categoric()

        self.children.update(
            select_genelists=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-genelists"),
                clearable=True, default=None, options=genelists,
                placeholder="Select Gene Lists", multi=True,
                style={"width": "100%"},
                options_callback=lambda: sorted(self.dataset.get_genelists()),
            ),
            select_clusterby=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-cluster_cells_by"),
                options=categoricals + ["barcode"], default=None, clearable=True,
                placeholder="Select Categorical (Optional)",
                options_callback=lambda: self.dataset.get_categoric() + ["barcode"],
            ),
            select_categoricals=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-categoricals"),
                options=categoricals, default=next(iter(categoricals), None),
                placeholder="Select Feature(s)", multi=True, clearable=True,
                style={"width": "100%"},
                options_callback=lambda: self.dataset.get_categoric(),
            ),
            select_layer=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-layer"),
                options=heatmap_params["layer"].allowed_values, default=heatmap_params["layer"].default,
            ),
            select_colormap=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-colormap"),
                options=heatmap_params["colormap"].allowed_values, default=heatmap_params["colormap"].default,
            ),
            select_genes=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-genes"),
                options=genes, default=None, clearable=True,
                placeholder="Select Genes", multi=True,
                style={"width": "100%"},
                options_callback=lambda: sorted(self.dataset.adata.var_names.tolist())
            ),
        )

        self.actions.update(
            plot_heatmap=PlotHeatmap(
                self.cid, dataset,
                select_genes_cid=self.children["select_genes"].cid,
                select_clusterby_cid=self.children["select_clusterby"].cid,
                select_categoricals_cid=self.children["select_categoricals"].cid,
                select_layer_cid=self.children["select_layer"].cid,
                select_colormap_cid=self.children["select_colormap"].cid,
            ),
            # edit_genelist=EditGeneList(self.cid, dataset),
            # add_genes_from_list=AddGenesFromList(CID(self.page_id, self.loc_class, "add_genes_from_list"), dataset),
        )

    def get_sidebar_params(self) -> list:
        divs = [
            html.Div([
                html.Label("Show Genes",className="param-label"),
                html.Div([
                    self.children["select_genes"].create_layout(),
                    self.children["select_genes"].get_stores(),
                    self.children["select_genelists"].create_layout(),
                    self.children["select_genelists"].get_stores(),
                ], style={"display": "flex", "gap": "10px"})
            ], className="param-row-stacked"),
            html.Div(
                children=[
                    html.Label("Cluster Group By", className="param-label"),
                    html.Div([
                        self.children["select_clusterby"].create_layout(),
                        self.children["select_clusterby"].get_stores(),
                    ], className="param-select"),
                ],
                className="param-row-stacked",
            ),
            html.Div(
                children=[
                    html.Label("Show Categorical Feature(s)", className="param-label"),
                    html.Div([
                        self.children["select_categoricals"].create_layout(),
                        self.children["select_categoricals"].get_stores(),
                    ], className="param-select"),
                ],
                className="param-row-stacked",
            ),
            html.Div(
                children=[
                    html.Label("Select Layer", className="param-label"),
                    html.Div([
                        self.children["select_layer"].create_layout(),
                        self.children["select_layer"].get_stores(),
                    ], className="param-select"),
                ],
                className="param-row-stacked",
            ),
            html.Div(
                children=[
                    html.Label("Select Colormap", className="param-label"),
                    html.Div([
                        self.children["select_colormap"].create_layout(),
                        self.children["select_colormap"].get_stores(),
                    ], className="param-select"),
                ],
                className="param-row-stacked",
            ),
        ]
        # divs.extend(components.params_layout(heatmap_params, f"{self.page_id}-heatmap"))


        return divs

    def create_layout(self):
        figure_layout = html.Div(
            children=[
                html.Div(
                    [
                        dcc.Loading(type="circle", children=[
                            html.Div([
                                dcc.Graph(
                                    id=self.cid.to_dict(), className=f"{self.loc_class}-plot"
                                )],
                            )]
                        )
                    ],
                    className=f"{self.loc_class}-body", id=f"heatmap-figure"
                )
            ]
        )

        return figure_layout
