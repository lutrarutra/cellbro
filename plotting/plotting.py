### Custom functions to plot SC RNA-seq data
import pandas as pd
import scanpy as sc
import numpy as np

import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rcParams

from sklearn.cluster import KMeans, AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram

def _dendrogram(model, **kwargs):
    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count        # axes[0][i+1].title.set_text(f"{clusterby.capitalize()} {groups[i]}")

    linkage_matrix = np.column_stack(
        [model.children_, np.linspace(0, 1, model.distances_.shape[0]+2)[1:-1], counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    return dendrogram(linkage_matrix, **kwargs)


def heatmap(adata, groupby, categorical_features, var_names=None, rank_genes_by=None, free_sort_cells=False, n_genes=10, sort_cells=True, sort_genes=True, quantiles=(0.0, 1.0), cmap="seismic", figsize=(20, None), dpi=50, fig_path=None):
    if isinstance(categorical_features, str):   
        categorical_features = [categorical_features]
    
    _grid = rcParams["axes.grid"]
    rcParams["axes.grid"] = False


    palettes = [sns.color_palette("tab10"), sns.color_palette("Paired"), sns.color_palette("Set2")]

    if var_names is None:
        rank_results = sc.tl.rank_genes_groups(
            adata, groupby=rank_genes_by if rank_genes_by else groupby,
            rankby_abs=True, method="t-test", copy=True
        ).uns["rank_genes_groups"]
        var_names = np.unique(np.array(list(map(list, zip(*rank_results["names"]))))[:,:n_genes].flatten())
    else:
        if isinstance(var_names, list):
            var_names = np.array(var_names)

    n_cat = len(categorical_features)
    h_cat = 0.5
    n_vars = len(var_names)
    h_vars = 0.3

    if figsize[1] is None:
        figsize = (figsize[0], n_cat * h_cat + n_vars * h_vars)

    r_cat = int(h_cat * 100.0 / (n_cat * h_cat + n_vars * h_vars))
    r_vars = 100 - r_cat

    f = plt.figure(figsize=figsize, dpi=dpi)

    gs = f.add_gridspec(
        1 + len(categorical_features), 2, hspace=0.2/figsize[1], wspace=0.01,
        height_ratios=[r_cat] * len(categorical_features) + [r_vars], width_ratios=[1, 20]
    )
    
    axes = gs.subplots(sharex="col", sharey="row")

    if sort_genes:
        model = AgglomerativeClustering(distance_threshold=0, n_clusters=None).fit(adata[:,var_names].X.T.toarray())
        gene_dendro = _dendrogram(model, ax=axes[-1, 0], orientation="right")
        gene_order = var_names[gene_dendro["leaves"]]

        icoord = np.array(gene_dendro["icoord"])
        icoord = icoord / (n_genes*10) * n_genes
        dcoord = np.array(gene_dendro["dcoord"])
        axes[-1, 0].clear()
        for xs, ys in zip(icoord, dcoord):
            axes[-1, 0].plot(ys, xs)
    else:
        gene_order = var_names

    if sort_cells:
        adata.obs["barcode"] = pd.Categorical(adata.obs.index)

        if free_sort_cells:
            if not f"dendrogram_barcode" in adata.uns.keys():
                sc.tl.dendrogram(adata, groupby="barcode", var_names=var_names)
            cell_dendro = adata.uns["dendrogram_barcode"]
            cell_order = cell_dendro["categories_ordered"]
        else:
            cell_order = []
            for cell_type in adata.obs[groupby].cat.categories.tolist():
                # Todo: when reculcustering, the order of the cells is not preserved
                if f"{cell_type}_order" in adata.uns.keys():
                    cell_order.extend(adata.uns[f"{cell_type}_order"])
                else:
                    dendro = sc.tl.dendrogram(adata[adata.obs[groupby] == cell_type], groupby="barcode", inplace=False)
                    cell_order.extend(dendro["categories_ordered"])
                    adata.uns[f"{cell_type}_order"] = dendro["categories_ordered"]

    else:
        cell_order = adata.obs.sort_values(groupby).index

    data = adata[cell_order, gene_order].layers["logcentered"].toarray().T
    vmin, vmax = np.quantile(data, q=quantiles)

    sns.heatmap(
        data, cmap=cmap, ax=axes[-1,-1],
        center=0, vmin=vmin, vmax=vmax, cbar=False, yticklabels=gene_order,
    )

    for i, categorical_feature in enumerate(categorical_features):
        palette = palettes[i % len(palettes)]
        samples = adata[cell_order, :].obs[categorical_feature].cat.codes
        clr = [palette[s % len(palette)] for s in samples]
        axes[i][1].vlines(np.arange(len(samples)), 0, 1, colors=clr, lw=5, zorder=10)
        axes[i][1].set_yticklabels([])
        axes[i][1].set_ylim([0,1])
        axes[i][1].patch.set_linewidth(2.0)
        axes[i][1].set_yticks([0.5])
        axes[i][1].set_yticklabels([categorical_feature.capitalize()])

        leg = f.legend(
            title=categorical_feature, labels=adata.obs[categorical_feature].cat.categories.tolist(),
            prop={"size": 24}, bbox_to_anchor=(0.95, 0.9 - 0.3*i), ncol=1, frameon=True, edgecolor="black",
            loc="upper left", facecolor="white"
        )

        plt.gca().add_artist(leg)
        palette = palettes[i%len(palettes)]
        for l, legobj in enumerate(leg.legendHandles):
            legobj.set_color(palette[l % len(palette)])
            legobj.set_linewidth(8.0)

    for ax in axes.flat:
        ax.yaxis.set_tick_params(length=0)
        ax.xaxis.set_tick_params(length=0)
        ax.set_xticklabels([])
    
    for ax in axes[:,0]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

    f.colorbar(
        plt.cm.ScalarMappable(norm=matplotlib.colors.TwoSlopeNorm(vmin=vmin, vmax=vmax, vcenter=0), cmap=cmap),
        ax=axes, orientation="vertical", fraction=0.05, pad=0.01, shrink=5.0 / figsize[1]
    )

    
    if fig_path:
        plt.savefig(fig_path, bbox_inches="tight")

    plt.show()
    rcParams["axes.grid"] = _grid