import multiprocessing

import scanpy as sc
import imgui

import analysis.plotting as pl

class Projection():
    def __init__(self, dataset, app):
        self.app = app
        self.type = "Projection"
        self.dataset = dataset
        self.finished = False
        self.leiden = True
        self.keys = sorted(list(dataset.adata.obs.columns)) + sorted(list(dataset.adata.var.index))
        self.selected_keys = []
        self.proposal_keys = list(dataset.adata.obs.columns)
        self.calc_params = {}
        self.plot_params = {
            "ncols": 1,
            "frameon": False,
        }

        self.leiden_params = {
            "resolution": 1.0,
            "random_state": 0,
        }

        self.query = ""

    def find_keys(self, query):
        keys = []
        for key in self.keys:
            if query.upper() in key.upper():
                keys.append(key)
                if len(keys) > 20:
                    return keys

        return keys

    def draw(self):
        if self.finished:
            return False, None
        
        imgui.begin(f"{self.type} Setup")
        
        imgui.text("Calculation Parameters")
        for key, value in self.calc_params.items():
            if isinstance(value, int):
                _, self.calc_params[key] = imgui.input_int(key, value)
            elif isinstance(value, float):
                _, self.calc_params[key] = imgui.input_float(key, value)

        imgui.dummy(0, 20)
        imgui.text("Plotting Parameters")

        imgui.text("Search:")
        query_changed, self.query  = imgui.input_text(
                "Key (Column/GeneName)", self.query, 256
        )
        if query_changed:
            self.proposal_keys = self.find_keys(self.query)

        imgui.begin_child("Available Keys text", 200, 20, border=False)
        imgui.text("Available Keys:")
        imgui.end_child()
        imgui.same_line()
        imgui.begin_child("Selected Keys text", 200, 20, border=False)
        imgui.text("Selected Keys:")
        imgui.end_child()

        imgui.begin_child("Available Keys", 200, -10, border=True)
        for key in self.proposal_keys:
            clicked, _ = imgui.selectable(key, key in self.selected_keys)
            if clicked:
                if key in self.selected_keys:
                    self.selected_keys.remove(key)
                else:
                    self.selected_keys.append(key)

        imgui.end_child()
        imgui.same_line()
        imgui.begin_child("Selected Keys", 200, -10, border=True)
        for key in self.selected_keys:
            if imgui.button(key):
                self.selected_keys.remove(key)

        imgui.end_child()

        for key, value in self.plot_params.items():
            if isinstance(value, int):
                _, self.plot_params[key] = imgui.input_int(key, value)
            elif isinstance(value, float):
                _, self.plot_params[key] = imgui.input_float(key, value)
            elif isinstance(value, bool):
                _, self.plot_params[key] = imgui.checkbox(key, value)

        imgui.dummy(0, 20)
        imgui.text("Leiden")
        imgui.same_line()
        _, self.leiden = imgui.checkbox("Enabled", self.leiden)
        if self.leiden:
            for key, value in self.leiden_params.items():
                if isinstance(value, int):
                    _, self.leiden_params[key] = imgui.input_int(key, value)
                elif isinstance(value, float):
                    _, self.leiden_params[key] = imgui.input_float(key, value)

        if imgui.button("Apply"):
            self.finished = True
            imgui.end()
            return False, self

        imgui.end()
        return True, None

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

        self.app.processes["umap_plot"] = multiprocessing.Process(target=pl.umap_projection, args=(self.app.dataset, self.plot_params,))
        self.app.processes["umap_plot"].start()