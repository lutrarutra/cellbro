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


class Projection():
    def __init__(self, dataset, params: ParamsDict):
        self.dataset = dataset
        self.params = params

    @classmethod
    def parse_params(cls, params):
        key = cls.get_key()
        print(key)
        projection_params = dict(
            [
                (param_key.replace(f"{key}_", ""), params[param_key])
                for param_key in params.keys()
                if param_key.startswith(f"{key}_")
            ]
        )
        print(projection_params)
        return projection_params

    @classmethod
    @abstractmethod
    def get_key(cls) -> str:
        ...

    @staticmethod
    @abstractmethod
    def get_type() -> ProjectionType:
        ...

    @staticmethod
    @abstractmethod
    def get_params() -> ParamsDict:
        ...

    @abstractmethod
    def apply(self) -> str:
        ...

    def add_params(self):
        if f"{self.get_type()}_params" not in self.dataset.adata.uns.keys():
            self.dataset.adata.uns[f"{self.get_type()}_params"] = {}

        rerun = (
            self.dataset.adata.uns[f"{self.get_type()}_params"] != self.params.unravel()
        )
        self.dataset.adata.uns[f"{self.get_type()}_params"] = self.params.unravel()

        return rerun

    @staticmethod
    def get_layout(projection_cls, page_id_prefix):
        divs = []
        for key, param in projection_cls.get_params().items():
            if param.type == bool:
                inp = dbc.Switch(
                    id=f"{page_id_prefix}-projection-{projection_cls.get_type().value}-{key}",
                    value=param.default,
                )
            else:
                inp = dbc.Input(
                    id=f"{page_id_prefix}-projection-{projection_cls.get_type().value}-{key}",
                    type=param.input_type,
                    value=param.value,
                    step=param.step if param.step != None else 0.1,
                    className="param-input",
                    placeholder=param.placeholder,
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


# def parse_params(params: dict):
#     projection_type = params.pop("projection_type")
#     projection_color = params.pop("projection_color")
#     key = None
#     if projection_type == "UMAP":
#         key = "umap"
#     elif projection_type == "t-SNE":
#         key = "tsne"
#     elif projection_type == "SCVI-UMAP":
#         key = "scvi_umap"
#     else:
#         key = "trimap"

#     projection_params = dict(
#         [
#             (param_key.replace(f"{key}_", ""), params[param_key])
#             for param_key in params.keys()
#             if param_key.startswith(f"{key}_")
#         ]
#     )
#     return projection_type, dict(color=projection_color, params=projection_params)
    # elif projection_type == "PCA":
    #     return PCA(dataset=dataset, color=color, params=params)
