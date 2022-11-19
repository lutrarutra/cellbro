import multiprocessing

import scanpy as sc

from .projection import Projection

class Trimap(Projection):
    def __init__(self, dataset, app):
        super().__init__(dataset, app)
        self.type = "Trimap"
        self.calc_params = dict(
            n_components=2, n_inliers=10, n_outliers=5, n_random=5
        )
        self.tl_documentation = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.tl.trimap.html"
        self.pl_documentation = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pl.trimap.html"

    def apply(self):
        sc.external.tl.trimap(
            self.dataset.adata, **self.calc_params
        )

        if self.leiden:
            self.selected_keys.append("leiden")
            sc.tl.leiden(
                self.dataset.adata, **self.leiden_params
            )
        
        self.plot_params["color"] = self.selected_keys

        self.app.processes["trimap_plot"] = multiprocessing.Process(target=sc.external.pl.trimap, args=(self.dataset.adata,), kwargs=self.plot_params)
        self.app.processes["trimap_plot"].start()