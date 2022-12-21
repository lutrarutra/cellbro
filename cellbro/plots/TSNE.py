import scanpy as sc

from cellbro.plots.Projection import Projection, ProjectionType
from cellbro.util.Param import Param, ParamsDict


class TSNE(Projection):
    def __init__(self, dataset, params):
        super().__init__(dataset, TSNE._params.update(params))

    def apply(self) -> str:
        sc.tl.tsne(self.dataset.adata, **self.params.unravel())
        return "X_tsne"

    @staticmethod
    def get_type() -> ProjectionType:
        return ProjectionType.TSNE

    @classmethod
    def get_key(cls) -> str:
        return "tsne"

    @staticmethod
    def get_params() -> ParamsDict:
        return TSNE._params

    _params = ParamsDict(
        [
            Param(
                key="n_pcs",
                name="Num. PCs",
                default=None,
                type=int,
                nullable=True,
                description="",
                placeholder="Optional",
                step=1,
            ),
            # Param(key="n_components", name="3D Projection", default=False, type=bool, nullable=False, description=""),
        ]
    )
