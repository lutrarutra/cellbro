import scanpy as sc
import pandas as pd

import util.FileFormat as FileFormat, util.FileType as FileType, util.IO as IO


class Dataset():
    def __init__(self, path, file_type, logger):
        self.path = path
        self.file_type = file_type
        self.preprocessed = False
        self.annotated = False

        if path.endswith(".csv"):
            self.file_format = FileFormat.CSV.value
        elif path.endswith(".tsv"):
            self.file_format = FileFormat.TSV.value
        elif path.endswith(".h5") or path.endswith(".hdf5"):
            self.file_format = FileFormat.H5.value
        elif path.endswith(".rds"):
            self.file_format = FileFormat.RDS.value
        else:
            logger.error(f"File format {'.' + path.split('.')[-1]} not supported...")
            return None

        self.logger = logger
        self.adata = IO.load_file(path, self.file_format, self.file_type, self.logger)
        # TODO: check duplicate barcodes
        self.adata.var_names_make_unique()
        self.adata.obs_names_make_unique()
        self.adata.var['mt'] = self.adata.var_names.str.startswith("MT.")
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    def annotate(self, path):
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

        
        