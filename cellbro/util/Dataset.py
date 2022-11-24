import scanpy as sc

class Dataset():
    def __init__(self, path):
        self.path = path
        print("Reading File")
        self.adata = sc.read_h5ad(self.path)
        print("Filtering Cells")
        sc.pp.filter_cells(self.adata, min_genes=3000)
        sc.pp.filter_genes(self.adata, min_cells=10)

        print("Calculating QC Metrics")
        self.adata.var["mt"] = self.adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=["mt"], percent_top=False, log1p=False, inplace=True)

        print("Normalizing Counts")
        self.adata.layers["counts"] = self.adata.X.copy()
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)

        print("Calculating Metrics")
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata, random_state=0)
        sc.tl.umap(self.adata, min_dist=0.3)
        print("Dataset Ready!")