import scanpy as sc

from cellbro.plots.Projection import Projection, ProjectionType
from cellbro.util.Param import Param, ParamsDict


class Trimap(Projection):
    def __init__(self, dataset, params):
        super().__init__(dataset, Trimap._params.update(params))

    def apply(self) -> str:
        params = self.params.unravel()
        if params.pop("3d_projection"):
            key = "X_trimap_3d"
            self.dataset.adata.obsm[key] = sc.external.tl.trimap(
                self.dataset.adata, copy=True, n_components=3, **params
            ).obsm["X_trimap"].copy()

        else:
            key = "X_trimap"
            sc.external.tl.trimap(self.dataset.adata, **params)

        return key

    @classmethod
    def get_key(cls) -> str:
        return "trimap"

    @staticmethod
    def get_type() -> ProjectionType:
        return ProjectionType.TRIMAP

    @staticmethod
    def get_params() -> ParamsDict:
        return Trimap._params

    _params = ParamsDict(
        [
            Param(
                key="n_inliers",
                name="Num. Inliers",
                default=10,
                type=int,
                description="",
                step=1,
            ),
            Param(
                key="n_outliers",
                name="Num. Outliers",
                default=5,
                type=int,
                description="",
                step=1,
            ),
            Param(
                key="n_random",
                name="Num. Random",
                default=5,
                type=int,
                description="",
                step=1,
            ),
            # Param(
            #     key="weight_adj",
            #     name="Weight Adj.",
            #     default=500.0,
            #     type=float,
            #     description="",
            #     step=50.0,
            # ),
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
