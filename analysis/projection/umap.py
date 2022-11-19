import multiprocessing

import scanpy as sc

from .projection import Projection

class UMAP(Projection):
    def __init__(self, dataset, app):
        super().__init__(dataset, app)
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
            self.dataset.adata,
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
            sc.tl.leiden(
                self.dataset.adata,
                resolution=self.leiden_params["resolution"],
                random_state=self.leiden_params["random_state"],
            )
        
        self.plot_params["color"] = self.selected_keys

        self.app.processes["umap_plot"] = multiprocessing.Process(target=sc.pl.umap, args=(self.dataset.adata,), kwargs=self.plot_params)
        self.app.processes["umap_plot"].start()