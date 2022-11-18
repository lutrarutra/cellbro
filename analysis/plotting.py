import scanpy as sc

import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rcParams
import seaborn as sns

def qc_mt_plot(dataset):
    sc.pl.scatter(dataset.adata, x="total_counts", y="pct_counts_mt", color="n_genes_by_counts")


def umap_projection(dataset, args):
    print(dataset.adata.shape)
    sc.pl.umap(dataset.adata, **args)