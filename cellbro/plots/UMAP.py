import plotly.express as px
import scanpy as sc

from cellbro.plots.Projection import Projection, ProjectionType
from cellbro.util.Param import Param, ParamsDict


class UMAP(Projection):
    def __init__(self, dataset, color, params):
        super().__init__(dataset, color, UMAP._params.update(params))

    def apply(self):
        sc.tl.umap(self.dataset.adata, **self.params.unravel())

    @staticmethod
    def get_type() -> ProjectionType:
        return ProjectionType.UMAP

    @staticmethod
    def get_key() -> str:
        return "X_umap"

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
                key="n_components",
                name="3D Projection",
                default=False,
                type=bool,
                nullable=False,
                description="",
            ),
        ]
    )