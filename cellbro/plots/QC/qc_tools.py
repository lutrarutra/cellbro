import plotly.graph_objects as go

from ...util.Param import Param, ParamsDict

import scanpy as sc
import scipy

default_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=5, b=5, l=5, r=5),
)

qc_params = ParamsDict(
    [
        Param(
            key="pct_counts_mt",
            name="Max MT%",
            default=5.0,
            type=float,
            description="",
            step=1,
        ),
        Param(
            key="min_genes",
            name="Min. Genes (per cell)",
            default=200,
            type=int,
            description="",
            step=1,
        ),
        Param(
            key="min_cells",
            name="Min. Cells (per gene)",
            default=3,
            type=int,
            description="",
            step=1,
        ),
    ]
)


def apply_dispersion_qc(dataset):
    if not "cv2" in dataset.adata.var.columns or not "mu" in dataset.adata.var.columns:
        ncounts = dataset.adata.layers["ncounts"]
        if isinstance(ncounts, scipy.sparse.csr_matrix):
            ncounts = ncounts.toarray()
        dataset.adata.var["cv2"] = (ncounts.std(0) / ncounts.mean(0)) ** 2
        dataset.adata.var["mu"] = ncounts.mean(0)


def apply_mt_qc(dataset):
    if not "pct_counts_mt" in dataset.adata.obs.columns:
        dataset.adata.var["mt"] = dataset.adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            dataset.adata, qc_vars=["mt"], percent_top=False, log1p=False, inplace=True
        )


def filter(self, submit):
    # Makes sure that filtering is not done on initial load
    if submit is None:
        return list(self.params.values())  # !=

    self.dataset.adata = self.dataset.adata[
        self.dataset.adata.obs.pct_counts_mt < self.params["pct_counts_mt"], :
    ].copy()
    sc.pp.filter_cells(self.dataset.adata, min_genes=self.params["min_genes"])
    sc.pp.filter_genes(self.dataset.adata, min_cells=self.params["min_cells"])

    return list(self.params.values())
