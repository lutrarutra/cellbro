import scanpy as sc

import imgui

import util.gui
class PreprocessProgress():
    def __init__(self):
        self.finished = False
        self.step = 0
        self.steps = ["Annotation", "Filtering", "Mitochondrial QC", "Normalization"]

    def draw(self):
        if self.finished:
            return False, None
        
        imgui.begin("Preprocessing Pipeline")
        for i, step in enumerate(self.steps):
            if i > self.step:
                imgui.push_style_color(imgui.COLOR_TEXT, 1.0, 0.0, 0.0)
            elif i == self.step:
                imgui.push_style_color(imgui.COLOR_TEXT, 1.0, 0.7, 0.0)
            else:
                imgui.push_style_color(imgui.COLOR_TEXT, 0.0, 1.0, 0.0)
            imgui.text(f"{i+1}. {step}")
            imgui.pop_style_color()
            
        imgui.end()
        return True, None

class AskPreprocessForm():
    def __init__(self):
        self.finished = False

    def draw(self):
        if self.finished:
            return False, None
        
        imgui.begin("Question")
        imgui.text("Do you want to preprocess the dataset?")
        if imgui.button("Yes"):
            self.finished = True
            imgui.end()
            return False, True
        imgui.same_line()
        if imgui.button("No"):
            self.finished = True
            imgui.end()
            return False, False

        imgui.end()
        return True, None

class AnnotateForm():
    def __init__(self):
        self.finished = False

    def draw(self):
        if self.finished:
            return False, None
        
        imgui.begin("Annotate")
        imgui.text("Do you have phenodata? (.csv/.tsv)")
        imgui.same_line()
        util.gui.tooltip("Comma/tab separated file with header as first row and cellbarcodes as first column.")
        if imgui.button("Yes"):
            imgui.end()
            return False, True
        imgui.same_line()
        if imgui.button("No"):
            imgui.end()
            return False, False
        imgui.end()
        return True, None

class FilterForm():
    def __init__(self, dataset):
        self.dataset = dataset
        self.finished = False
        self.min_genes = 200
        self.min_counts = 0
        self.min_cells = 3

    def draw(self):
        if self.finished:
            return False, None
        
        imgui.begin("Preprocessing: Filter Genes and Cells")
        # sc.pp.filter_cells(...)
        imgui.text(f"Filter Cells: {self.dataset.adata.shape[0]}")
        min_genes_changed, self.min_genes = imgui.input_int("Min genes", self.min_genes)
        util.gui.tooltip("Minimum number of genes expressed required for a cell to pass filtering.")
        if min_genes_changed:
            self.min_counts = 0

        imgui.text("or")
        min_counts_changed, self.min_counts = imgui.input_int("Min counts", self.min_counts)
        util.gui.tooltip("Maximum number of counts required for a cell to pass filtering.")
        if min_counts_changed:
            self.min_genes = 0

        imgui.dummy(0, 20)

        # sc.pp.filter_genes(...)
        imgui.text(f"Filter Genes: {self.dataset.adata.shape[1]}")
        _, self.min_cells = imgui.input_int("Min. cells", self.min_cells)
        util.gui.tooltip("Minimum number of cells expressed required for a gene to pass filtering.")
        
        if imgui.button("Apply"):
            if self.min_genes > 0:
                sc.pp.filter_cells(self.dataset.adata, min_genes=self.min_genes)
            elif self.min_counts > 0:
                sc.pp.filter_cells(self.dataset.adata, min_counts=self.min_counts)
            sc.pp.filter_genes(self.dataset.adata, min_cells=self.min_cells)
            self.finished = True
            imgui.end()
            return False, None
        imgui.same_line()
        if imgui.button("Skip"):
            self.finished = True
            imgui.end()
            return False, None

        imgui.end()
        return True, None


class MTQCForm():
    def __init__(self, dataset):
        self.dataset = dataset
        self.finished = False
        self.pct_counts_mt = 5

    def draw(self):
        if self.finished:
            return False, None
        
        imgui.begin("Preprocessing: Remove cells with high mitochondrial content")

        _, self.pct_counts_mt = imgui.input_float("Percentage", self.pct_counts_mt)
        util.gui.tooltip("Max percentage of mitochondrial gene counts.")

        if imgui.button("Apply"):
            self.dataset.adata = self.dataset.adata[self.dataset.adata.obs.pct_counts_mt < self.pct_counts_mt, :]
            self.finished = True
            imgui.end()
            return False, None

        imgui.end()
        return True, None

class NormalizeForm():
    def __init__(self, dataset):
        self.dataset = dataset
        self.finished = False
        self.target_sum = -1
        self.target_sum_mean = True

    def draw(self):
        if self.finished:
            return False, None
        
        imgui.begin("Preprocessing: Normalize")
        check_box_clicked, self.target_sum_mean = imgui.checkbox("Mean normalization", self.target_sum_mean)
        if check_box_clicked:
            self.target_sum = -1
        util.gui.tooltip("If None, after normalization, each observation (cell) has a total count equal\nto the median of total counts for observations (cells) before normalization.")
        
        imgui.text("or")

        target_sum_changed, self.target_sum = imgui.input_int("Total counts", self.target_sum)
        if target_sum_changed:
            self.target_sum_mean = False

        util.gui.tooltip("Total counts per cell.")

        if imgui.button("Apply"):
            self.dataset.adata.layers["counts"] = self.dataset.adata.X.copy()
            sc.pp.normalize_total(self.dataset.adata, target_sum=None if self.target_sum_mean else self.target_sum)

            self.finished = True
            imgui.end()
            return False, None

        imgui.end()
        return True, None

