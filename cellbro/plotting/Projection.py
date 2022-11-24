
import plotly.graph_objects as go
import plotly.express as px
from dash import html, dcc, Input, Output, State

import scanpy as sc

from cellbro.util.Param import *

projection_layout = go.Layout(
    paper_bgcolor='white',
    plot_bgcolor='white',
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
)

umap_params = ParamsDict([
    Param(key="min_dist", name="Min. Distance", default=1.0, type=float, description="", step=0.1, _min=0.0),
    Param(key="spread", name="Spread", default=1.0, type=float, description="", step=0.1, _min=0.0),
])

tsne_params = ParamsDict([
    Param(key="n_pcs", name="Num. PCs", default=None, type=int, nullable=True, description=""),
])

trimap_params = ParamsDict([
    Param(key="n_inliers", name="Num. Inliers", default=10, type=int, description="", step=1),
    Param(key="n_outliers", name="Num. Outliers", default=5, type=int, description="", step=1),
    Param(key="n_random", name="Num. Random", default=5, type=int, description="", step=1),
    Param(key="weight_adj", name="Weight Adj.", default=500.0, type=float, description="", step=50.0),
])

pca_params = ParamsDict([])

class Projection():
    def __init__(self, key, dataset, color):
        self.key = key
        self.dataset = dataset
        self.color = color

    # Outputs from _plot()
    @staticmethod
    def get_callback_outputs():
        return [
            Output(component_id="projection-plot", component_property="figure"),
            Output(component_id="projection-umap", component_property="style"),
            Output(component_id="projection-tsne", component_property="style"),
            Output(component_id="projection-trimap", component_property="style"),
            Output(component_id="projection-pca", component_property="style")
        ]
        
    # Inputs to Projection
    @staticmethod
    def get_callback_inputs():
        inputs = {
            "projection_color": Input(component_id="projection-color", component_property="value"),
            "projection_type": Input(component_id="projection-type", component_property="value")
        }
        for key in umap_params.keys():
            inputs[f"umap_{key}"] = Input(component_id=f"projection-umap-{key}", component_property="value")

        for key in tsne_params.keys():
            inputs[f"tsne_{key}"] = Input(component_id=f"projection-tsne-{key}", component_property="value")
        
        for key in trimap_params.keys():
            inputs[f"trimap_{key}"] = Input(component_id=f"projection-trimap-{key}", component_property="value")
        
        for key in pca_params.keys():
            inputs[f"pca_{key}"] = Input(component_id=f"projection-pca-{key}", component_property="value")

        return inputs
        

    @staticmethod
    def create_layout(dataset):
        layout = html.Div([
            html.Div(
                children=[
                    html.Div(children=[
                    # Projection type celect
                        html.Div(children=[
                            html.Label("Projection Type"),
                            dcc.Dropdown(["UMAP", "Trimap", "t-SNE", "PCA"], value="UMAP", id="projection-type", clearable=False),
                        ], style={"padding": "10px", "flex": "1"}),

                        # Projection Hue celect
                        html.Div(children=[
                            html.Label("Color"),
                            dcc.Dropdown(dataset.adata.obs_keys(), value=dataset.adata.obs_keys()[0], id="projection-color", clearable=False),
                        ], style={"padding": "10px", "flex": "1"}),
                    ], style={"padding": "10px", "flex": "1", "display": "flex", "flex-direction": "row"}),
                    html.Div(id="projection-parameters", children=[dcc.Loading(
                            id="loading-projection", type="circle",
                            children=[Projection.params_layout()],
                    )]),
                ], style={"padding": "10px", "flex": "1"}
            ),
            html.Div(children=[
                dcc.Loading(
                    id="loading-projection", type="circle",
                    children=[html.Div(dcc.Graph(id="projection-plot"))],
            )], style={"padding": "10px", "flex": "3"}),

        ], style={"display": "flex", "flex-direction": "row"})
        return layout

    @staticmethod
    def _param_layout(projection_type, params):
        layout = html.Div(children=[
            html.Div(children=[
                html.Label(
                    param.name,
                    style={"flex": "1"}
                ),
                dcc.Input(
                    id=f"projection-{projection_type}-{key}", type=param.input_type, value=param.value,
                    step=param.step if param.step != None else 0.1, style={"flex": "1"}
                ),
            ], style={"display":"flex", "padding": "10px", "flex": "1"}) for key, param in params.items()
        ], style={"display": "none"}, id=f"projection-{projection_type.lower()}")
        return layout

    @staticmethod
    def params_layout():
        return html.Div(children=[
            Projection._param_layout("umap", umap_params),
            Projection._param_layout("tsne", tsne_params),
            Projection._param_layout("trimap", trimap_params),
            Projection._param_layout("pca", pca_params),
        ])

    def plot(self):
        projection_figure = px.scatter(
            x=self.dataset.adata.obsm[self.key][:,0], y=self.dataset.adata.obsm[self.key][:,1],
            color=self.dataset.adata.obs[self.color], color_discrete_sequence=sc.pl.palettes.default_20,
            labels={"x": f"{self.type} 1", "y": f"{self.type} 2"}
        )
        projection_figure.update_layout(projection_layout)
        projection_figure.update_xaxes(
            scaleanchor="y", scaleratio=1,
        )
        return projection_figure

    def add_params(self, adata):
        if f"{self.type}_params" not in adata.uns.keys():
            adata.uns[f"{self.type}_params"] = {}

        rerun = (adata.uns[f"{self.type}_params"] != self.params.unravel())
        if rerun:
            print(type(adata.uns[f"{self.type}_params"]), type(self.params))
            print(self.params.unravel())
            print(adata.uns[f"{self.type}_params"]) 
        adata.uns[f"{self.type}_params"] = self.params.unravel()

        return rerun

class UMAP(Projection):
    def __init__(self, dataset, color, params):
        super().__init__("X_umap", dataset, color)
        self.type = "UMAP"
        self.params = umap_params.update(params)
        rerun = self.add_params(dataset.adata)
        if self.key not in dataset.adata.obsm.keys() or rerun:
            sc.tl.umap(dataset.adata, **self.params.unravel())

        self.add_params(dataset.adata)

class TSNE(Projection):
    def __init__(self, dataset, color, params):
        super().__init__("X_tsne", dataset, color)
        self.type = "t-SNE"
        self.params = tsne_params.update(params)
        rerun = self.add_params(dataset.adata)
        if self.key not in dataset.adata.obsm.keys() or rerun:
            sc.tl.tsne(dataset.adata, **self.params.unravel())

        self.add_params(dataset.adata)

class PCA(Projection):
    def __init__(self, dataset, color, params):
        super().__init__("X_pca", dataset, color)
        self.type = "PCA"
        self.params = pca_params.update(params)
        rerun = self.add_params(dataset.adata)
        if self.key not in dataset.adata.obsm.keys() or rerun:
            sc.tl.pca(dataset.adata, **self.params.unravel())

        self.add_params(dataset.adata)


class Trimap(Projection):
    def __init__(self, dataset, color, params):
        super().__init__("X_trimap", dataset, color)
        self.type = "Trimap"
        self.params = trimap_params.update(params)
        rerun = self.add_params(dataset.adata)
        if self.key not in dataset.adata.obsm.keys() or rerun:
            sc.external.tl.trimap(dataset.adata, **self.params.unravel())
        
        

