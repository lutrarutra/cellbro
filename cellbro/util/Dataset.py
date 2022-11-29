import scanpy as sc

import scout

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
        scout.tl.scale_log_center(self.adata, target_sum=None)

        ncounts = self.adata.layers["ncounts"].toarray()
        self.adata.var["cv2"] = (ncounts.std(0)/ncounts.mean(0))**2
        self.adata.var["mu"] = ncounts.mean(0)

        print("Calculating Metrics")
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata, random_state=0)


        print("Dataset Ready!")