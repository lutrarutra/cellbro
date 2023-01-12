# cellbrowser

## Quality Control
![alt text](https://github.com/lutrarutra/scout/blob/main/figures/qc.png?raw=true)

## Cell Browser with Projections, Violins, and Heatmap
![alt text](https://github.com/lutrarutra/scout/blob/main/figures/cells.png?raw=true)

## Differential Expression (DE)
![alt text](https://github.com/lutrarutra/scout/blob/main/figures/de.png?raw=true)
![alt text](https://github.com/lutrarutra/scout/blob/main/figures/de_settings.png?raw=true)

## Gene Set Enrichment Analysis (GSEA)
![alt text](https://github.com/lutrarutra/scout/blob/main/figures/gsea.png?raw=true)
![alt text](https://github.com/lutrarutra/scout/blob/main/figures/gsea_settings.png?raw=true)

## PCA with Projection, Correlation Circle, Correlation w.r.t. cell features, and Explained Variance
![alt text](https://github.com/lutrarutra/scout/blob/main/figures/pca.png?raw=true)

### How to Run?
1. Download source:
    - `git clone --recurse-submodules https://github.com/lutrarutra/cellbro`
    - `cd cellbro`
2. Install required libraries:
    - `pip install -r requirements.txt`
3. Run
    - `python cellbro/run.py <path_to_file.h5ad>
