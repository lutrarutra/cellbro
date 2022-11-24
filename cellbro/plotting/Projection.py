
import plotly.graph_objects as go
import plotly.express as px

import scanpy as sc

projection_layout = go.Layout(
    paper_bgcolor='white',
    plot_bgcolor='white',
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
)

class Projection():
    def __init__(self, key, dataset, color):
        self.key = key
        self.dataset = dataset
        self.color = color

    def plot(self):
        projection_figure = px.scatter(
            x=self.dataset.adata.obsm[self.key][:,0], y=self.dataset.adata.obsm[self.key][:,1],
            color=self.dataset.adata.obs[self.color], color_discrete_sequence=sc.pl.palettes.default_20,
            labels={"x": f"{self.type} 1", "y": f"{self.type} 2"}
        )
        projection_figure.update_layout(projection_layout)
        projection_figure.update_xaxes(
            scaleanchor="y", scaleratio=1,
        )
        return projection_figure

class UMAP(Projection):
    def __init__(self, dataset, color):
        super().__init__("X_umap", dataset, color)
        self.type = "UMAP"
        if self.key not in dataset.adata.obsm.keys():
            sc.tl.umap(dataset.adata)

class TSNE(Projection):
    def __init__(self, dataset, color):
        super().__init__("X_tsne", dataset, color)
        self.type = "t-SNE"
        if self.key not in dataset.adata.obsm.keys():
            sc.tl.tsne(dataset.adata)

class PCA(Projection):
    def __init__(self, dataset, color):
        super().__init__("X_pca", dataset, color)
        self.type = "PCA"
        if self.key not in dataset.adata.obsm.keys():
            sc.tl.pca(dataset.adata)
    
class Trimap(Projection):
    def __init__(self, dataset, color):
        super().__init__("X_trimap", dataset, color)
        self.type = "Trimap"
        if self.key not in dataset.adata.obsm.keys():
            sc.external.tl.trimap(dataset.adata)

