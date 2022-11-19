import multiprocessing

import scanpy as sc
import imgui

from .projection import Projection

class PCA(Projection):
    def __init__(self, dataset, app):
        super().__init__(dataset, app)
        self.type = "PCA"
        self.calc_params = dict()
        self.dimensions = [(0,1)]
        self.current_tab = 1
        self.tl_documentation = "https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.pca.html"
        self.pl_documentation = "https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.pca.html"

    def ask_dimensions(self):
        if imgui.button("Add dimension"):
            self.dimensions.append((0,1))

        imgui.begin_child("Available Keys text", imgui.get_window_width()-10, 200, border=False)
        imgui.push_item_width((imgui.get_window_width()-150)*0.5)

        for i, (x,y) in enumerate(self.dimensions):
            changed, val = imgui.input_int(f"[{i+1}] X", x)
            if changed and val >= 0:
                self.dimensions[i] = (val, y)
            imgui.same_line()
            changed, val = imgui.input_int(f"[{i+1}] Y", y)
            if changed and val >= 0:
                self.dimensions[i] = (x, val)

        imgui.pop_item_width()
        imgui.end_child()

    def apply(self):
        if self.leiden:
            self.selected_keys.append("leiden")
            sc.tl.leiden(
                self.dataset.adata, **self.leiden_params
            )
        
        self.plot_params["color"] = self.selected_keys
        self.plot_params["dimensions"] = list(set(self.dimensions))

        self.app.processes["pca_plot"] = multiprocessing.Process(target=sc.pl.pca, args=(self.dataset.adata,), kwargs=self.plot_params)
        self.app.processes["pca_plot"].start()