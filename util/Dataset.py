import os, datetime

import imgui

import scanpy as sc
import pandas as pd

import util.FileFormat as FileFormat, util.FileType as FileType, util.IO as IO


class Dataset():
    def __init__(self, path, file_type, logger):
        self.path = path
        self.file_type = file_type
        self.preprocessed = False
        self.annotated = False
        self.annotation_files = []

        if path.endswith(".csv"):
            self.file_format = FileFormat.CSV
        elif path.endswith(".tsv"):
            self.file_format = FileFormat.TSV
        elif path.endswith(".h5") or path.endswith(".hdf5") or path.endswith(".h5ad"):
            self.file_format = FileFormat.H5
        elif path.endswith(".rds"):
            self.file_format = FileFormat.RDS
        else:
            logger.error(f"File format {'.' + path.split('.')[-1]} not supported...")
            return None

        self.logger = logger

    def load_file(self):
        self.adata = IO.load_file(self.path, self.file_format, self.file_type, self.logger)
        # TODO: check duplicate barcodes
        self.adata.var_names_make_unique()
        self.adata.obs_names_make_unique()
        self.adata.var['mt'] = self.adata.var_names.str.startswith("MT-") | self.adata.var_names.str.startswith("MT.")
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    def draw(self):
        layers = list(self.adata.layers.keys())
        n_features = max(len(layers), max(len(self.adata.obs.columns), len(self.adata.var.columns)))

        imgui.text(f"Cells: {self.adata.shape[0]}")
        imgui.same_line()
        imgui.dummy(10, 0)
        imgui.same_line()
        imgui.text(f"Genes: {self.adata.shape[1]}")

        imgui.dummy(0, 20)

        imgui.text("Dataset Content")
        imgui.columns(3, "contentlist")
        imgui.text("Cell Features")
        imgui.next_column()
        imgui.text("Gene Features")
        imgui.next_column()
        imgui.text("Layers")
        imgui.separator()
        imgui.next_column()
        imgui.push_style_color(imgui.COLOR_TEXT, 0.7, 0.7, 0.7, 1.0)
        for i in range(n_features):
            if i < len(self.adata.obs.columns):
                imgui.text(f"{self.adata.obs.columns[i]} ({type(self.adata.obs.iloc[0,i]).__name__})")
            else:
                imgui.text("")
            imgui.next_column()
            if i < len(self.adata.var.columns):
                imgui.text(f"{self.adata.var.columns[i]} ({type(self.adata.var.iloc[0,i]).__name__})")
            else:
                imgui.text("")
            imgui.next_column()
            if i < len(self.adata.layers.keys()):
                imgui.text(layers[i])
            else:
                imgui.text("")
            imgui.separator()
            imgui.next_column()

        imgui.pop_style_color()
        imgui.columns(1)

        imgui.dummy(0, 20)

        imgui.text("Loaded Files:")
        imgui.columns(4, "filelist")

        imgui.text("Type")
        imgui.next_column()
        imgui.text("File")
        imgui.next_column()
        imgui.text("Size")
        imgui.next_column()
        imgui.text("Last Modified")
        imgui.next_column()

        imgui.separator()
        
        imgui.text(f"SC Data")
        imgui.next_column()
        imgui.text(os.path.basename(self.path))
        imgui.next_column()
        imgui.text(f"{os.stat(self.path).st_size / (1024 * 1024):.1f} MB")
        imgui.next_column()
        imgui.text(datetime.datetime.fromtimestamp(os.path.getmtime(self.path)).strftime("%Y-%m-%d %H:%M:%S"))

        imgui.separator()
        imgui.push_style_color(imgui.COLOR_TEXT, 0.7, 0.7, 0.7, 1.0)
        for annotation in self.annotation_files:
            imgui.next_column()
            imgui.text(f"Annotation")
            imgui.next_column()
            imgui.text(os.path.basename(annotation))
            imgui.next_column()
            imgui.text(f"{os.stat(annotation).st_size / (1024 * 1024):.1f} MB")
            imgui.next_column()
            imgui.text(datetime.datetime.fromtimestamp(os.path.getmtime(annotation)).strftime("%Y-%m-%d %H:%M:%S"))
            imgui.separator()

        imgui.pop_style_color()
        imgui.columns(1)

    def annotate(self, path):
        self.annotation_files.append(path)
        sep = "," if path.split(".")[-1] == "csv" else "\t"
        annotation = pd.read_csv(path, sep=sep, index_col=0)

        if annotation.shape[0] != self.adata.shape[0]:
            self.logger.warning(f"Phenodata file does not contain all barcodes ({annotation.shape[0]} != {self.adata.shape[0]}), absent barcodes are filled with NaNs.")

        if set(annotation.index) != set(self.adata.obs.index):
            self.logger.warning(f"Phenodata file does not contain all barcodes: ({len(list(set(annotation.index)).symmetric_difference(set(self.adata.obs.index)))}/{self.adata.obs.index.shape[0]}), absent barcodes are filled with NaNs.")

        for col in annotation.columns:
            self.adata.obs[col] = None
            self.adata.obs.loc[annotation.index, col] = annotation[col]
            self.logger.info(f"Added column {col} to adata.obs")

    def post_process_init(self):
        self.adata.layers["ncounts"] = self.adata.X.copy()
        sc.pp.log1p(self.adata)
        self.adata.layers["centered"] = self.adata.layers["ncounts"] - self.adata.layers["ncounts"].mean(axis=0)
        self.adata.layers["logcentered"] = self.adata.X - self.adata.X.mean(axis=0)
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata, n_neighbors=15, random_state=0)
        self.preprocessed = True


        
        