from dataclasses import dataclass

import cellbro.plots.QC as QC

import scanpy as sc
import pandas as pd
import scipy

import scout

@dataclass
class Dataset:
    path: str
    adata: sc.AnnData

    def __init__(self, path):
        self.path = path
        print("Reading File")

        self.adata = sc.read_h5ad(self.path)

        self.adata.obs["barcode"] = pd.Categorical(self.adata.obs_names)

        if "gene_lists" not in self.adata.uns.keys():
            self.adata.uns["gene_lists"] = {}

        self.adata.uns["scvi_setup_params"] = {}
        self.adata.uns["scvi_model_params"] = {}
        self.adata.uns["scvi_train_params"] = {}

        if "log1p" in self.adata.uns.keys():
            self.adata.uns["log1p"] = {"base": None}

        QC.apply_dispersion_qc(self)
        QC.apply_mt_qc(self)
        
        print("Dataset Ready!")


    def pp(self):
        print("Filtering Cells")
        sc.pp.filter_cells(self.adata, min_genes=3000)
        sc.pp.filter_genes(self.adata, min_cells=10)

        print("Normalizing Counts")
        scout.tl.scale_log_center(self.adata, target_sum=None)

        print("Calculating Metrics")
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata, random_state=0)

    def qc_done(self):
        if not "pct_counts_mt" in self.adata.obs.columns:
            return False
        
        if not "cv2" in self.adata.var.columns:
            return False
        
        if not "mu" in self.adata.var.columns:
            return False

        return True

    def get_categoric(self):
        return list(
            set(self.adata.obs.columns)
            - set(self.adata.obs._get_numeric_data().columns) - set(["barcode"])
        )

    def get_obs_features(self, include_gene_lists=False):
        res = sorted(list(self.adata.obs_keys()))
        
        if include_gene_lists:
            for gene_list in self.get_gene_lists():
                res.append(f"Gene List: {gene_list}")

        res.extend(sorted(list(self.adata.var_names)))

        if "barcode" in res:
            res.remove("barcode")


        return res

    def get_numeric(self):
        return list(self.adata.obs._get_numeric_data().columns)

    def get_neighbors(self):
        return [key for key in self.adata.uns_keys() if "neighbors" in key]

    def get_gene_lists(self, gene=None):
        if "gene_lists" not in self.adata.uns.keys():
            return []
        if gene is None:
            return list(self.adata.uns["gene_lists"].keys())

        return [
            gl
            for gl in self.adata.uns["gene_lists"].keys()
            if gene in self.adata.uns["gene_lists"][gl]
        ]

    def get_genes_from_list(self, gene_list):
        return self.adata.uns["gene_lists"][gene_list]

    def update_gene_lists(self, gene, gene_lists):
        for gl_key in self.get_gene_lists():
            if gl_key in gene_lists:
                if gene not in self.adata.uns["gene_lists"][gl_key]:
                    self.adata.uns["gene_lists"][gl_key].append(gene)
            else:
                if gene in self.adata.uns["gene_lists"][gl_key]:
                    self.adata.uns["gene_lists"][gl_key].remove(gene)

        return self.get_gene_lists(gene)

    def get_rank_genes_groups(self):
        return sorted([key[11:] for key in self.adata.uns_keys() if "rank_genes_" in key])

    def get_available_refs(self, groupby):
        return sorted(list(self.adata.uns["rank_genes_" + groupby].keys()))

    def get_scvi_projections(self):
        return [x for x in list(self.adata.obsm.keys()) if "X_umap_scvi" in x]

    def add_gsea_result(self, res, groupby, reference):
        if "gsea" not in self.adata.uns.keys():
            self.adata.uns["gsea"] = {}

        if groupby not in self.adata.uns["gsea"].keys():
            self.adata.uns["gsea"][groupby] = {}

        self.adata.uns["gsea"][groupby][reference] = res
