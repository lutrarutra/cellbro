import multiprocessing, webbrowser, threading

import pandas as pd
import scanpy as sc
import numpy as np

import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rcParams

from plotting import Figure

import time

class Label():
    def __init__(self, ax, dx=20):
        self.annotation = ax.annotate(
            "", xy=(0,0), xytext=(20,20), textcoords="offset points",
            bbox=dict(boxstyle="round", fc="w"),
            arrowprops=dict(arrowstyle="->")
        )
        self.annotation.get_bbox_patch().set_facecolor("whitesmoke")
        self.annotation.get_bbox_patch().set_alpha(0.8)

    def set_visible(self, visible):
        self.annotation.set_visible(visible)

    def remove(self):
        self.annotation.remove()


class LivePlot(Figure.Figure):
    def __init__(self, app, type, plot_fnc):
        super().__init__(app, type, self.plot)
        self.fig = None
        self.ax = None
        self.path_collection = None
        self.annotation = None
        self.labels = {}
        self.plot_fnc = plot_fnc
        self.selected_label_key = None
        # Not working properly when using "repeat keys"-accessbility setting in OS
        self.key_down = False

    def _update(self, ind, label):
        pos = self.path_collection.get_offsets()[ind["ind"][0]]
        label.annotation.xy = pos
        text = "\n".join([
            f"{self.annotation.index[n]}: (logFC: {self.annotation['logFC'][n]:.1f}, pval: {self.annotation['-log_pvals_adj'][n]:.1f})" for n in ind["ind"]
        ])
        label.annotation.set_text(text)

    def _add_label(self, ind, key):
        self.labels[key] = Label(
            self.ax
        )
        self.labels[key].set_visible(True)
        self._update(ind, self.labels[key])
        self.fig.canvas.draw_idle()

    def _remove_label(self, key, draw=True):
        self.labels[key].set_visible(False)
        del self.labels[key]
        if draw:
            self.fig.canvas.draw_idle()

    def _move_label(self, key, dx, dy):
        x, y = self.labels[key].annotation.get_position()
        self.labels[key].annotation.set_position((x-dx, y-dy))
        self.fig.canvas.draw_idle()

    def on_key_pressed(self, event):
        if self.selected_label_key is None:
            return

        if self.key_down:
            return

        if event.inaxes == self.ax:
            _, ind = self.path_collection.contains(event)
        else:
            return        

        if event.key == "left":
            self._move_label(self.selected_label_key, 5, 0)
        elif event.key == "right":
            self._move_label(self.selected_label_key, -5, 0)
        elif event.key == "up":
            self._move_label(self.selected_label_key, 0, -5)
        elif event.key == "down":
            self._move_label(self.selected_label_key, 0, 5)
        else:
            return
            
        self.key_down = True
        

    def on_key_released(self, event):
        if event.key in ["left", "right", "up", "down"]:
            self.key_down = False


    def on_mouse_pressed(self, event):
        if (event.button == matplotlib.backend_bases.MouseButton.LEFT):
            if event.inaxes == self.ax:
                cont, ind = self.path_collection.contains(event)
            else:
                return

            if not cont:
                return

            key = ind["ind"].tobytes()
            if key not in self.labels.keys():
                self._add_label(ind, key)
                self.selected_label_key = key
            else:
                self._remove_label(key)
                self.selected_label_key = None

        elif (event.button == matplotlib.backend_bases.MouseButton.RIGHT):
            if event.inaxes == self.ax:
                cont, ind = self.path_collection.contains(event)
                if cont:
                    webbrowser.open(f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={self.annotation.index[ind['ind'][0]]}")


    def _plot(self, params):
        self.fig, self.ax, self.path_collection, self.annotation = self.plot_fnc(**params)
        self.fig.canvas.mpl_connect("button_press_event", self.on_mouse_pressed)
        self.fig.canvas.mpl_connect("key_press_event", self.on_key_pressed)
        self.fig.canvas.mpl_connect("key_release_event", self.on_key_released)
        plt.show()

    def plot(self, params):
        self.app.figures[self.id] = multiprocessing.Process(target=self._plot, args=(params,))
        self.app.figures[self.id].start()
        