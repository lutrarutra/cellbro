
from cellbro.plots.Projection import Projection, ProjectionType
from cellbro.util.Param import *

import scanpy as sc

class MDE(Projection):
    def __init__(self, dataset, color, params):
        super().__init__(dataset, color, MDE._params.update(params))

    def apply(self):
        sc.tl.umap(self.dataset.adata, **self.params.unravel())

    @staticmethod
    def get_type() -> ProjectionType:
        return ProjectionType.UMAP

    @staticmethod
    def get_key() -> str:
        return "X_mde"

    @staticmethod
    def get_params() -> ParamsDict:
        return MDE._params

    _params = ParamsDict([])
