import webbrowser, multiprocessing

import imgui

import scanpy as sc
import pandas as pd

from util import Query


class Violin():
    def __init__(self, app):
        self.app = app
        self.params = dict(
            keys=[],
            groupby=None,
            inner="box",
            stripplot=True,
            jitter=0.4,
            rotation=0,
        )

        self.inner_proposal = ["None", "box", "quartile", "point", "stick"]
        self.inner_selected = 1

        self.key_query = Query.Query(
            sorted(list(self.app.dataset.adata.obs.columns)) + sorted(list(self.app.dataset.adata.var.index)),
            proposal_keys=list(self.app.dataset.adata.obs.columns)
        )
        self.selected_keys = []

        self.proposal_groupby = ["None"] + [
            x for x in self.app.dataset.adata.obs.columns \
                if type(self.app.dataset.adata.obs.dtypes[x]) == pd.CategoricalDtype \
                    or type(self.app.dataset.adata.obs.dtypes[x]) == str
        ]
        self.selected_groupby = "None"
        self.query = ""

    def draw(self):
        _, self.params["stripplot"] = imgui.checkbox("Stripplot", self.params["stripplot"])
        imgui.push_item_width(imgui.get_window_width()*0.3)
        if self.params["stripplot"]:
            imgui.same_line()
            _, self.params["jitter"] = imgui.input_float("Jitter", self.params["jitter"])
        
        # TODO: inner not working or removed?
        # clicked, self.inner_selected = imgui.combo("Inner Plot", self.inner_selected, self.inner_proposal)
        _, self.params["rotation"] = imgui.input_int("Label Rotation", self.params["rotation"], step=45)
        imgui.pop_item_width()

        imgui.text("Group by:")
        imgui.begin_child("search groupby", imgui.get_window_width()-20, 200, border=True)
        for key in self.proposal_groupby:
            if imgui.radio_button(key, key == self.selected_groupby):
                self.selected_groupby = key
        imgui.end_child()

        imgui.dummy(0, 20)

        imgui.text("Search:")
        imgui.push_item_width(imgui.get_window_width()*0.5)
        query_changed, self.query  = imgui.input_text(
                "Key (Feature/GeneName)", self.query, 256
        )
        imgui.pop_item_width()
        if query_changed:
            self.key_query.query(self.query)

        imgui.begin_child("Available features text", (imgui.get_window_width()-30)*0.5, 40, border=False)
        imgui.text("Available Features:")
        imgui.end_child()
        imgui.same_line()
        imgui.begin_child("Selected features text", (imgui.get_window_width()-30)*0.5, 40, border=False)
        imgui.text("Selected Features:")
        imgui.end_child()
        imgui.begin_child("Available Features", (imgui.get_window_width()-30)*0.5, -50, border=True)
        for key in self.key_query.proposal_keys:
            clicked, _ = imgui.selectable(key, key in self.selected_keys)
            if clicked:
                if key in self.selected_keys:
                    self.selected_keys.remove(key)
                else:
                    self.selected_keys.append(key)

        imgui.end_child()
        imgui.same_line()
        imgui.begin_child("Selected Features", (imgui.get_window_width()-30)*0.5, -50, border=True)
        for key in self.selected_keys:
            if imgui.button(key):
                self.selected_keys.remove(key)

        imgui.end_child()

        imgui.set_cursor_pos((imgui.get_cursor_pos()[0], imgui.get_window_height() - 40))
        if imgui.button("Plot"):
            self.apply()
            return False
        imgui.same_line()
        if imgui.button("Cancel"):
            return False
        imgui.same_line()
        if imgui.button("Documentation"):
            webbrowser.open("https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.violin.html")

        return True

    def apply(self):
        self.params["inner"] = self.inner_proposal[self.inner_selected] if not "None" else None
        self.params["keys"] = self.selected_keys
        self.params["groupby"] = self.selected_groupby if self.selected_groupby != "None" else None

        self.params["adata"] = self.app.dataset.adata
        self.plot(self.params)

        # self.app.figures["violin_plot"] = multiprocessing.Process(target=sc.pl.violin, args=(self.app.dataset.adata,), kwargs=self.params)
        # self.app.figures["violin_plot"].start()
