import plotly.express as px
import scanpy as sc

from cellbro.plots.Projection import Projection, ProjectionType
from cellbro.util.Param import Param, ParamsDict


class UMAP(Projection):
    def __init__(self, dataset, params: ParamsDict):
        super().__init__(dataset, UMAP._params.update(params))

    def apply(self) -> str:
        params = self.params.unravel()
        if params.pop("3d_projection"):
            key = "X_umap_3d"
            self.dataset.adata.obsm[key] = sc.tl.umap(
                self.dataset.adata, copy=True, n_components=3, **params
            ).obsm["X_umap"].copy()
        else:
            key = "X_umap"
            sc.tl.umap(self.dataset.adata, **params)
        
        return key

    @classmethod
    def get_key(cls) -> str:
        return "umap"

    @staticmethod
    def get_type() -> ProjectionType:
        return ProjectionType.UMAP

    @staticmethod
    def get_params() -> ParamsDict:
        return UMAP._params

    _params = ParamsDict(
        [
            Param(
                key="min_dist",
                name="Min. Distance",
                default=1.0,
                type=float,
                description="",
                step=0.1,
                _min=0.0,
            ),
            Param(
                key="spread",
                name="Spread",
                default=1.0,
                type=float,
                description="",
                step=0.1,
                _min=0.0,
            ),
            Param(
                key="3d_projection",
                name="3D Projection",
                default=False,
                type=bool,
                nullable=False,
                description="",
            ),
        ]
    )


class SCVI_UMAP(UMAP):
    def __init__(self, dataset, params: ParamsDict):
        super().__init__(dataset, params)

    def apply(self):
        params = self.params.unravel()
        if params.pop("3d_projection"):
            key = "X_umap_scvi_3d"

            self.dataset.adata.obsm[key] = sc.tl.umap(
                self.dataset.adata, neighbors_key="neighbors_scvi", n_components=3,
                copy=True, **params
            ).obsm["X_umap"].copy()

            return key
        else:
            key = "X_umap_scvi"

            self.dataset.adata.obsm[key] = sc.tl.umap(
                self.dataset.adata, neighbors_key="neighbors_scvi",
                copy=True, **params
            ).obsm["X_umap"].copy()

            return key

    @classmethod
    def get_key(cls) -> str:
        return "scvi_umap"

    @staticmethod
    def get_type() -> ProjectionType:
        return ProjectionType.SCVI_UMAP
