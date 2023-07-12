# cellbrowser

## Quality Control
![alt text](https://github.com/lutrarutra/cellbro/blob/main/figures/qc.png?raw=true)

## Cell Browser with Projections, Violins, and Heatmap
![alt text](https://github.com/lutrarutra/cellbro/blob/main/figures/cells.png?raw=true)
![alt text](https://github.com/lutrarutra/cellbro/blob/main/figures/projection.png?raw=true)
![alt text](https://github.com/lutrarutra/cellbro/blob/main/figures/heatmap.png?raw=true)

## Differential Expression (DE)
![alt text](https://github.com/lutrarutra/cellbro/blob/main/figures/de.png?raw=true)
![alt text](https://github.com/lutrarutra/cellbro/blob/main/figures/de_settings.png?raw=true)

## Gene Set Enrichment Analysis (GSEA)
![alt text](https://github.com/lutrarutra/cellbro/blob/main/figures/gsea.png?raw=true)
![alt text](https://github.com/lutrarutra/cellbro/blob/main/figures/gsea_settings.png?raw=true)

## PCA with Projection, Correlation Circle, Correlation w.r.t. cell features, and Explained Variance
![alt text](https://github.com/lutrarutra/cellbro/blob/main/figures/pca.png?raw=true)

### How to Run?
1. Download source:
    - `git clone https://github.com/lutrarutra/cellbro`
    - `git clone https://github.com/lutrarutra/scout`
2. Install required libraries:
    - `pip install -r cellbro/requirements.txt`
    - `pip install -e ./scout`
    - pip install -e ./cellbro`
3. Run (h5ad file expects specific format: `preprocess_example.py`)
    - `python cellbro/run.py <path_to_file.h5ad>
