import multiprocessing

import scanpy as sc

from .Projection import Projection
from plotting import Figure
class TSNE(Projection, Figure.Figure):
    def __init__(self, app):
        Projection.__init__(self, app)
        Figure.Figure.__init__(self, app, "tsne", sc.pl.tsne)
        self.type = "t-SNE"
        self.calc_params = dict(
            perplexity=30, early_exaggeration=12
        )
        self.tl_documentation = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.tsne.html"
        self.pl_documentation = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.tsne.html"


    def apply(self):
        sc.tl.tsne(
            self.app.dataset.adata, **self.calc_params
        )

        if self.leiden:
            self.selected_keys.append("leiden")
            self._leiden.apply()

        self.plot_params["color"] = self.selected_keys
        self.plot_params["adata"] = self.app.dataset.adata

        self.plot(self.plot_params)