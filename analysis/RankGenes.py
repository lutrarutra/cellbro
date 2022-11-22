import imgui

import scanpy as sc
import pandas as pd

from util import Task

class RankGenes():
    def __init__(self, app):
        self.app = app
        self.params = dict(
            groupby=None,
            method=None,
            corr_method=None,
            reference=None,
        )
        self.available_groupby = [
            x for x in self.app.dataset.adata.obs.columns \
                if type(self.app.dataset.adata.obs.dtypes[x]) == pd.CategoricalDtype \
                    or type(self.app.dataset.adata.obs.dtypes[x]) == str
        ]
        self.i_groupby = 0

        self.available_methods = ["t-test", "logreg", "wilcoxon", "t-test_overestim_var"]
        self.i_method = 0

        self.available_corr_methods = ["benjamini-hochberg", "bonferroni"]
        self.i_corr_method = 0

        self.available_reference = ["rest"] + self.available_groupby
        self.i_reference = 0

        self.documentation = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html"

    def draw(self):
        _, self.i_groupby = imgui.combo("Groupby", self.i_groupby, self.available_groupby)
        _, self.i_method = imgui.combo("Method", self.i_method, self.available_methods)
        _, self.i_corr_method = imgui.combo("Correction Method", self.i_corr_method, self.available_corr_methods)
        _, self.i_reference = imgui.combo("Reference", self.i_reference, self.available_reference)

    def apply(self, start_thread=False):
        self.params["groupby"] = self.available_groupby[self.i_groupby]
        self.params["method"] = self.available_methods[self.i_method]
        self.params["corr_method"] = self.available_corr_methods[self.i_corr_method]
        self.params["reference"] = self.available_reference[self.i_reference]
        print(self.app.dataset.adata.uns.keys())
        if start_thread:
            self.app.task_handler.add_task(
                "rank_genes_groups",
                Task.Task(target=sc.tl.rank_genes_groups, args=(self.app.dataset.adata,), kwargs=self.params)
            )
            self.app.task_handler.tasks["rank_genes_groups"].start()
        else:
            sc.tl.rank_genes_groups(self.app.dataset.adata, **self.params)