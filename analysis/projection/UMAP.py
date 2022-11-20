import multiprocessing

import scanpy as sc

from .Projection import Projection
from analysis import Figure
class UMAP(Projection, Figure.Figure):
    def __init__(self, app):
        Projection.__init__(self, app)
        Figure.Figure.__init__(self, app, "umap", sc.pl.umap)
        self.type = "UMAP"
        self.calc_params = {
            "min_dist" : 0.5,
            "spread" : 1.0,
            "n_components" : 2,
            "alpha" : 1.0,
            "gamma" : 1.0,
            "negative_sample_rate" : 5,
            "random_state" : 0,
        }
        self.tl_documentation = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.umap.html"
        self.pl_documentation = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.umap.html"

    def apply(self):
        sc.tl.umap(
            self.app.dataset.adata,
            min_dist=self.calc_params["min_dist"],
            spread=self.calc_params["spread"],
            n_components=self.calc_params["n_components"],
            alpha=self.calc_params["alpha"],
            gamma=self.calc_params["gamma"],
            negative_sample_rate=self.calc_params["negative_sample_rate"],
            random_state=self.calc_params["random_state"],
        )
        if self.leiden:
            self.selected_keys.append("leiden")
            self._leiden.apply()
        
        self.plot_params["color"] = self.selected_keys
        self.plot_params["adata"] = self.app.dataset.adata

        self.plot(self.plot_params)