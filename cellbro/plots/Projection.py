from abc import ABC, abstractmethod
from enum import Enum

import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import Input, Output, State, dcc, html

from cellbro.util.Param import *

projection_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
    margin=dict(t=10, b=10, l=10, r=10),
)


# pca_params = ParamsDict(
#     [
#         Param(
#             key="n_components",
#             name="3D Projection",
#             default=False,
#             type=bool,
#             nullable=False,
#             description="",
#         ),
#     ]
# )


class ProjectionType(Enum):
    UMAP = "umap"
    TRIMAP = "trimap"
    TSNE = "tsne"
    MDE = "mde"
    SCVI_UMAP = "scvi_umap"
    # SCVI_MDE = "scvi_mde"
    # PCA = "PCA"


class Projection(ABC):
    def __init__(self, dataset, color, params: ParamsDict):
        self.dataset = dataset
        self.color_label = color
        self.color = color

        # TODO: something smarter :)
        if "n_components" in params.params.keys():
            params.params["n_components"].value = 3 if params.params["n_components"].value else 2

        self.params = params

        if color in self.dataset.adata.obs_keys():
            self.color = self.dataset.adata.obs[color]
        else:
            self.color = (
                self.dataset.adata.X[:, self.dataset.adata.var.index.get_loc(color)]
                .toarray()
                .T[0]
            )

        rerun = self.add_params()
        if self.get_key() not in self.dataset.adata.obsm.keys() or rerun:
            self.apply()

    @staticmethod
    @abstractmethod
    def get_type() -> ProjectionType:
        ...

    @staticmethod
    @abstractmethod
    def get_key() -> str:
        ...

    @staticmethod
    @abstractmethod
    def get_params() -> ParamsDict:
        ...

    @abstractmethod
    def apply(self):
        ...

    def plot(self):
        if (
            "n_components" not in self.params.keys()
            or self.params["n_components"].value == 2
        ):
            fig = px.scatter(
                x=self.dataset.adata.obsm[self.get_key()][:, 0],
                y=self.dataset.adata.obsm[self.get_key()][:, 1],
                color=self.color,
                color_discrete_sequence=sc.pl.palettes.default_20,
                labels={
                    "x": f"{self.get_type().value} 1",
                    "y": f"{self.get_type().value} 2",
                },
            )
        else:
            fig = px.scatter_3d(
                x=self.dataset.adata.obsm[self.get_key()][:, 0],
                y=self.dataset.adata.obsm[self.get_key()][:, 1],
                z=self.dataset.adata.obsm[self.get_key()][:, 2],
                color=self.color,
                color_discrete_sequence=sc.pl.palettes.default_20,
                labels={
                    "x": f"{self.get_type().value} 1",
                    "y": f"{self.get_type().value} 2",
                    "z": f"{self.get_type().value} 3",
                },
            )

            fig.update_traces(marker=dict(size=2))

        projection_layout["legend"] = dict(title=self.color_label.capitalize())
        fig.update_layout(projection_layout)
        fig.update_xaxes(
            scaleanchor="y",
            scaleratio=1,
        )
        return fig

    def add_params(self):
        if f"{self.get_type()}_params" not in self.dataset.adata.uns.keys():
            self.dataset.adata.uns[f"{self.get_type()}_params"] = {}

        rerun = (
            self.dataset.adata.uns[f"{self.get_type()}_params"] != self.params.unravel()
        )
        self.dataset.adata.uns[f"{self.get_type()}_params"] = self.params.unravel()

        return rerun

    @staticmethod
    def get_layout(projection_cls):
        divs = []
        for key, param in projection_cls.get_params().items():
            if param.type == bool:
                inp = dbc.Switch(
                    id=f"projection-{projection_cls.get_type().value}-{key}",
                    value=param.default,
                )
            else:
                inp = dcc.Input(
                    id=f"projection-{projection_cls.get_type().value}-{key}",
                    type=param.input_type,
                    value=param.value,
                    step=param.step if param.step != None else 0.1,
                    className="param-input",
                )

            divs.append(
                html.Div(
                    children=[
                        html.Label(
                            param.name,
                            className="param-label",
                        ),
                        inp,
                    ],
                    className="param-row",
                )
            )


        return divs
# class PCA(Projection):
#     def __init__(self, dataset, color, params):
#         super().__init__("X_pca", dataset, color)
#         self.type = "PCA"
#         params["n_components"] = 3 if params["n_components"] else 2
#         self.params = pca_params.update(params)
#         rerun = self.add_params(dataset.adata)
#         # if self.key not in dataset.adata.obsm.keys() or rerun:
#         #     sc.tl.pca(dataset.adata, **self.params.unravel())

#         self.add_params(dataset.adata)


def parse_params(params):
    projection_type = params.pop("projection_type")
    projection_color = params.pop("projection_color")
    key = None
    if projection_type == "UMAP":
        key = "umap"
    elif projection_type == "t-SNE":
        key = "tsne"
    elif projection_type == "SCVI-UMAP":
        key = "scvi_umap"
    else:
        key = "trimap"

    projection_params = dict(
        [
            (param_key.replace(f"{key}_", ""), params[param_key])
            for param_key in params.keys()
            if param_key.startswith(f"{key}_")
        ]
    )
    print(projection_params)
    return projection_type, dict(color=projection_color, params=projection_params)
    # elif projection_type == "PCA":
    #     return PCA(dataset=dataset, color=color, params=params)
