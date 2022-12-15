from dataclasses import dataclass

import scanpy as sc

import scout


@dataclass
class Dataset:
    path: str
    adata: sc.AnnData

    def __init__(self, path):
        self.path = path
        print("Reading File")
        self.adata = sc.read_h5ad(self.path)
        print("Filtering Cells")
        sc.pp.filter_cells(self.adata, min_genes=3000)
        sc.pp.filter_genes(self.adata, min_cells=10)

        print("Calculating QC Metrics")
        self.adata.var["mt"] = self.adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            self.adata, qc_vars=["mt"], percent_top=False, log1p=False, inplace=True
        )

        print("Normalizing Counts")
        scout.tl.scale_log_center(self.adata, target_sum=None)

        ncounts = self.adata.layers["ncounts"].toarray()
        self.adata.var["cv2"] = (ncounts.std(0) / ncounts.mean(0)) ** 2
        self.adata.var["mu"] = ncounts.mean(0)

        print("Calculating Metrics")
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata, random_state=0)
        self.adata.uns["gene_lists"] = {}
        print("Dataset Ready!")

        self.adata.uns["scvi_setup_params"] = {}
        self.adata.uns["scvi_model_params"] = {}
        self.adata.uns["scvi_train_params"] = {}

    def get_categoricals(self):
        return list(
            set(self.adata.obs.columns)
            - set(self.adata.obs._get_numeric_data().columns)
        )

    def get_continuous(self):
        return list(self.adata.obs._get_numeric_data().columns)

    def get_neighbors(self):
        return [key for key in self.adata.uns_keys() if "neighbors" in key]

    def get_gene_lists(self, gene=None):
        if "gene_lists" not in self.adata.uns:
            return []
        if gene is None:
            return list(self.adata.uns["gene_lists"].keys())

        return [
            gl
            for gl in self.adata.uns["gene_lists"].keys()
            if gene in self.adata.uns["gene_lists"][gl]
        ]

    def update_gene_lists(self, gene, gene_lists):
        for gl_key in self.get_gene_lists():
            if gl_key in gene_lists:
                if gene not in self.adata.uns["gene_lists"][gl_key]:
                    self.adata.uns["gene_lists"][gl_key].append(gene)
            else:
                if gene in self.adata.uns["gene_lists"][gl_key]:
                    self.adata.uns["gene_lists"][gl_key].remove(gene)

        return self.get_gene_lists(gene)
