import multiprocessing

import scanpy as sc

from .projection import Projection
class TSNE(Projection):
    def __init__(self, app):
        super().__init__(app)
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
            sc.tl.leiden(
                self.app.dataset.adata, **self.leiden_params
            )
        
        self.plot_params["color"] = self.selected_keys

        self.app.processes["tsne_plot"] = multiprocessing.Process(target=sc.pl.tsne, args=(self.app.dataset.adata,), kwargs=self.plot_params)
        self.app.processes["tsne_plot"].start()