import scanpy as sc

import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rcParams
import seaborn as sns

def qc_mt_plot(dataset):
    sc.pl.scatter(dataset.adata, x="total_counts", y="pct_counts_mt", color="n_genes_by_counts")


def umap_projection(dataset, args):
    sc.pl.umap(dataset.adata, **args)

def trimap_projection(dataset, args):
    sc.external.pl.trimap(dataset.adata, **args)

def tsne_projection(dataset, args):
    sc.pl.tsne(dataset.adata, **args)

def pca_projection(dataset, args):
    sc.pl.pca(dataset.adata, **args)