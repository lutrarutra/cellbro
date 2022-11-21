import imgui

import scanpy as sc
import pandas as pd

from plotting import Figure

from plotting import plotting as pl


class Heatmap(Figure.Figure):
    def __init__(self, app):
        super().__init__(app, "heatmap", pl.heatmap)
        self.app = app
        self.plot_params = dict(
            groupby=None,
            categorical_features=[],
            var_names=None,
            rank_genes_by=None,
            free_sort_cells=False,
            n_genes=10,
            sort_cells=True,
            sort_genes=True,
            quantiles=(0.0, 1.0),
        )

        self.i_groupby = 0

        self.selected_categoricals = []
        self.available_categoricals = [
            x for x in self.app.dataset.adata.obs.columns \
                if type(self.app.dataset.adata.obs.dtypes[x]) == pd.CategoricalDtype \
                    or type(self.app.dataset.adata.obs.dtypes[x]) == str
        ]

    def draw(self):
        _, self.i_groupby = imgui.combo("Groupby", self.i_groupby, self.available_categoricals)

        imgui.begin_child("Available categoricals text", (imgui.get_window_width()-30)*0.5, 40, border=False)
        imgui.text("Available Categorical Features:")
        imgui.end_child()
        imgui.same_line()
        imgui.begin_child("Selected categoricals text", (imgui.get_window_width()-30)*0.5, 40, border=False)
        imgui.text("Selected Categorical Features:")
        imgui.end_child()
        imgui.begin_child("Available categoricals Features", (imgui.get_window_width()-30)*0.5, 300, border=True)
        for key in self.available_categoricals:
            clicked, _ = imgui.selectable(key, key in self.selected_categoricals)
            if clicked:
                if key in self.selected_categoricals:
                    self.selected_categoricals.remove(key)
                else:
                    self.selected_categoricals.append(key)

        imgui.end_child()
        imgui.same_line()
        imgui.begin_child("Selected categoricals Features", (imgui.get_window_width()-30)*0.5, 300, border=True)
        for key in self.selected_categoricals:
            if imgui.button(key):
                self.selected_categoricals.remove(key)
        imgui.end_child()

        for key, value in self.plot_params.items():
            if isinstance(value, bool):
                _, self.plot_params[key] = imgui.checkbox(key, value)
            elif isinstance(value, int):
                _, self.plot_params[key] = imgui.input_int(key, value)
            elif isinstance(value, float):
                _, self.plot_params[key] = imgui.input_float(key, value)

        imgui.set_cursor_pos((imgui.get_cursor_pos()[0], imgui.get_window_height() - 40))
        if len(self.selected_categoricals) == 0:
            imgui.push_style_color(imgui.COLOR_TEXT, 0.5, 0.5, 0.5)
        else:
            imgui.push_style_color(imgui.COLOR_TEXT, 1, 1, 1)
        if imgui.button("Plot"):
            if len(self.selected_categoricals) > 0:
                self.apply()
                imgui.pop_style_color()
                return True
        imgui.pop_style_color()
        imgui.same_line()
        if imgui.button("Cancel"):
            return False

        return True

    def apply(self):
        self.plot_params["categorical_features"] = self.selected_categoricals
        self.plot_params["groupby"] = self.available_categoricals[self.i_groupby]
        self.plot_params["adata"] = self.app.dataset.adata
        self.plot(self.plot_params)
