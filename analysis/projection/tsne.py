import multiprocessing

import scanpy as sc

from .projection import Projection
from analysis import Figure
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
            sc.tl.leiden(
                self.app.dataset.adata, **self.leiden_params
            )
        
        self.plot_params["color"] = self.selected_keys
        self.plot_params["adata"] = self.app.dataset.adata
        self.plot(self.plot_params)
        # self.app.figures["tsne_plot"] = multiprocessing.Process(target=sc.pl.tsne, args=(self.app.dataset.adata,), kwargs=self.plot_params)
        # self.app.figures["tsne_plot"].start()