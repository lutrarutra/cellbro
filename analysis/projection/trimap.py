import multiprocessing
import random
import string

import scanpy as sc

from .projection import Projection
from analysis import Figure
class Trimap(Projection, Figure.Figure):
    def __init__(self, app):
        Projection.__init__(self, app)
        Figure.Figure.__init__(self, app, "trimap", sc.external.pl.trimap)
        self.type = "Trimap"
        self.calc_params = dict(
            n_components=2, n_inliers=10, n_outliers=5, n_random=5
        )
        self.tl_documentation = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.tl.trimap.html"
        self.pl_documentation = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pl.trimap.html"

    def apply(self):
        sc.external.tl.trimap(
            self.app.dataset.adata, **self.calc_params
        )

        if self.leiden:
            self.selected_keys.append("leiden")
            sc.tl.leiden(
                self.app.dataset.adata, **self.leiden_params
            )
        
        self.plot_params["color"] = self.selected_keys
        self.plot_params["adata"] = self.app.dataset.adata
        self.plot(self.plot_params)
        
        # self.app.figures[figure_id] = multiprocessing.Process(target=sc.external.pl.trimap, args=(self.app.dataset.adata,), kwargs=self.plot_params)
        # self.app.figures[figure_id].start()