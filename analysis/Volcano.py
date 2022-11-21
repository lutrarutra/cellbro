import webbrowser

import imgui

import scanpy as sc
import pandas as pd

from util import Query
from plotting import Figure
from plotting import LivePlot
from analysis import RankGenes

from plotting import plotting as pl
from util import Task


class Volcano(LivePlot.LivePlot):
    def __init__(self, app):
        super().__init__(app, "volcano", pl.volcano_plot)
        self.app = app
        self.rank_genes = RankGenes.RankGenes(app)
        self.plot_params = dict(
            hue=None,
            significance_threshold=0.05,
            vmin=0.0,
            vmax=0.0,
            vcenter=0.0,
        )

        self.auto_vmin = True
        self.auto_vmax = True
        self.auto_vcenter = True

        self.current_tab = 0
        
        self.available_hue = ["log_mu_expression", "mu_expression"] + list(self.app.dataset.adata.var.columns)
        self.i_hue = 0

    def apply(self):
        self.plot_params["hue"] = self.available_hue[self.i_hue]
        self.plot_params["adata"] = self.app.dataset.adata
        if self.auto_vcenter:
            self.plot_params["vcenter"] = None
        if self.auto_vmin:
            self.plot_params["vmin"] = None
        if self.auto_vmax:
            self.plot_params["vmax"] = None
        self.plot(self.plot_params)

    def draw(self):
        imgui.columns(2, "tabs")
        if imgui.selectable("Rank Genes", self.current_tab == 0)[0]:
            self.current_tab = 0
        imgui.next_column()
        if not "rank_genes_groups" in self.app.dataset.adata.uns.keys():
            imgui.push_style_color(imgui.COLOR_TEXT, 0.5, 0.5, 0.5)
            imgui.text("Plotting")
            imgui.pop_style_color()
        else:
            if imgui.selectable("Plotting", self.current_tab == 1)[0]:
                self.current_tab = 1
        imgui.columns(1)
        imgui.separator()

        if self.current_tab == 0:
            self.rank_genes.draw()
            imgui.set_cursor_pos((imgui.get_cursor_pos()[0], imgui.get_window_height() - 40))
            if imgui.button("Apply"):
                self.rank_genes.apply(start_thread=True)
                self.current_tab = 1
        elif self.current_tab == 1:
            _, self.i_hue = imgui.combo("Hue", self.i_hue, self.available_hue)
            _, self.plot_params["significance_threshold"] = imgui.input_float("Significance Threshold", self.plot_params["significance_threshold"])
            

            _, self.auto_vmin = imgui.checkbox("Auto", self.auto_vmin)
            if not self.auto_vmin:
                imgui.same_line()
                _, self.plot_params["vmin"] = imgui.input_float("vmin", self.plot_params["vmin"])

            _, self.auto_vcenter = imgui.checkbox("Auto", self.auto_vcenter)
            if not self.auto_vcenter:
                imgui.same_line()
                _, self.plot_params["vcenter"] = imgui.input_float("vcenter", self.plot_params["vcenter"])

            _, self.auto_vmax = imgui.checkbox("Auto", self.auto_vmax)
            if not self.auto_vmax:
                imgui.same_line()   
                _, self.plot_params["vmax"] = imgui.input_float("vmax", self.plot_params["vmax"])

            imgui.set_cursor_pos((imgui.get_cursor_pos()[0], imgui.get_window_height() - 40))
            if imgui.button("Plot"):
                self.apply()
                return True

        imgui.same_line()
        if imgui.button("Cancel"):
            return False

        if self.current_tab == 0:
            imgui.same_line()
            if imgui.button("Documentation"):
                webbrowser.open("https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html")

        return True