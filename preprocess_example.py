import tqdm

import scanpy as sc
import scout

# Data from 
# http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

adata = sc.read_10x_mtx("data/filtered_gene_bc_matrices/hg19/", var_names="gene_symbols", cache=True)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)

scout.tl.scale_log_center(adata, target_sum=None)

sc.pp.pca(adata)
sc.pp.neighbors(adata, random_state=0)
sc.tl.leiden(adata, random_state=0)
sc.tl.umap(adata)

cats = scout.tl.get_categoric(adata)
for cat in cats:
    scout.tl.rank_marker_genes(adata, groupby=cat)

adata.uns["gsea"] = {}
for cat in cats:
    adata.uns["gsea"][cat] = {}
    pbar = tqdm.tqdm(adata.uns["de"][cat].items())
    for key, de_df in pbar:
        adata.uns["gsea"][cat][key] = {}
        for geneset in ["KEGG_2021_Human", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021"]:
            pbar.set_postfix({"Category": cat, "target": key, "geneset": geneset})
            adata.uns["gsea"][cat][key][geneset] = scout.tl.GSEA(de_df, gene_set=geneset, lead_genes_type="str")

adata.write_h5ad("data/adata.h5ad")