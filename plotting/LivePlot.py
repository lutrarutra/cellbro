import multiprocessing, webbrowser

import pandas as pd
import scanpy as sc
import numpy as np

import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rcParams

from plotting import Figure


class Label():
    def __init__(self, ax, persistent=False, xytext=None):
        xytext = (20,20) if xytext is None else xytext
        self.annotation = ax.annotate(
            "", xy=(0,0), xytext=xytext, textcoords="offset points",
            bbox=dict(boxstyle="round", fc="w"),
            arrowprops=dict(arrowstyle="->")
        )
        self.persistent = persistent

    def set_visible(self, visible):
        self.annotation.set_visible(visible)

    def remove(self):
        self.annotation.remove()


class LivePlot(Figure.Figure):
    def __init__(self, app, type, plot_fnc):
        super().__init__(app, type, self.plot)
        self.fig = None
        self.ax = None
        self.plt = None
        self.annotation = None
        self.labels = {}
        self.plot_fnc = plot_fnc

    def update(self, ind, label):
        pos = self.plt.get_offsets()[ind["ind"][0]]
        label.xy = pos
        text = "\n".join([
            f"{self.annotation.index[n]}: (logFC: {self.annotation['logFC'][n]:.1f}, pval: {self.annotation['-log_pvals_adj'][n]:.1f})" for n in ind["ind"]
        ])
        label.set_text(text)
        label.get_bbox_patch().set_alpha(0.4)

    def on_click(self, event):
        if event.button == matplotlib.backend_bases.MouseButton.LEFT:
            if event.inaxes == self.ax:
                cont, ind = self.plt.contains(event)
            else:
                return

            if not cont:
                return

            key = ind["ind"].tobytes()
            if key not in self.labels.keys():
                xytext = (20,20) if self.annotation['logFC'][ind['ind'][0]] > 0 else (-120,20)
                self.labels[key] = Label(self.ax, persistent=True, xytext=xytext)
                self.fig.canvas.draw_idle()
                return

            self.labels[key].persistent = not self.labels[key].persistent

            self.fig.canvas.draw_idle()

        elif event.button == matplotlib.backend_bases.MouseButton.RIGHT:
            if event.inaxes == self.ax:
                cont, ind = self.plt.contains(event)
            else:
                return

            if not cont:
                return

            webbrowser.open(f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={self.annotation.index[ind['ind'][0]]}")

    def on_hover(self, event):
        if event.inaxes == self.ax:
            cont, ind = self.plt.contains(event)
        else:
            return

        key = ind["ind"].tobytes()
        if cont:
            if not key in self.labels.keys():
                xytext = (20,20) if self.annotation['logFC'][ind['ind'][0]] > 0 else (-120,20)
                self.labels[key] = Label(self.ax, persistent=False, xytext=xytext)
                self.update(ind, self.labels[key].annotation)
                self.labels[key].set_visible(True)
        else:
            for key in list(self.labels.keys()):
                if not self.labels[key].persistent:
                    self.labels[key].remove()
                    del self.labels[key]

        self.fig.canvas.draw_idle()

    def _plot(self, params):
        self.fig, self.ax, self.plt, self.annotation = self.plot_fnc(**params)
        self.fig.canvas.mpl_connect("motion_notify_event", self.on_hover)
        self.fig.canvas.mpl_connect("button_press_event", self.on_click)
        plt.show()

    def plot(self, params):
        self.app.figures[self.id] = multiprocessing.Process(target=self._plot, args=(params,))
        self.app.figures[self.id].start()
        