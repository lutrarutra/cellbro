import scanpy as sc

from cellbro.plots.Projection import Projection, ProjectionType
from cellbro.util.Param import Param, ParamsDict


class Trimap(Projection):
    def __init__(self, dataset, color, params):
        super().__init__(dataset, color, Trimap._params.update(params))

    def apply(self):
        sc.external.tl.trimap(self.dataset.adata, **self.params.unravel())

    @staticmethod
    def get_type() -> ProjectionType:
        return ProjectionType.TRIMAP

    @staticmethod
    def get_key() -> str:
        return "X_trimap"

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
            Param(
                key="weight_adj",
                name="Weight Adj.",
                default=500.0,
                type=float,
                description="",
                step=50.0,
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
