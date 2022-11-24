import tkinter as tk
from tkinter import filedialog

from dash import Dash, html, dcc, Input, Output, State
import plotly.express as px
import plotly.graph_objects as go

from cellbro.util.Dataset import Dataset
from cellbro.plotting.Projection import Projection, UMAP, TSNE, PCA, Trimap

import scanpy as sc

app = Dash(__name__)
app.enable_dev_tools(debug=True, dev_tools_hot_reload=False)

# root = tk.Tk()
# root.withdraw()
# path = filedialog.askopenfilename()
dataset = Dataset("data/vas.h5ad")

app.layout = html.Div([
    Projection.create_layout(dataset)
])

@app.callback(
    output=Projection.get_callback_outputs(),
    inputs=Projection.get_callback_inputs(),
)
def update_projection(projection_color, projection_type, **kwargs):
    
    if projection_type == "UMAP":
        projection_params = dict(
            [(key.replace("umap_", ""), kwargs[key]) for key in kwargs.keys() if key.startswith("umap_")]
        )
        return UMAP(dataset, projection_color, projection_params).plot(), {"display": "block"}, {"display": "none"}, {"display": "none"}, {"display": "none"}

    elif projection_type == "t-SNE":
        projection_params = dict(
            [(key.replace("tsne_", ""), kwargs[key]) for key in kwargs.keys() if key.startswith("tsne_")]
        )
        return TSNE(dataset, projection_color, projection_params).plot(), {"display": "none"}, {"display": "block"}, {"display": "none"}, {"display": "none"}
    
    elif projection_type == "Trimap":
        projection_params = dict(
            [(key.replace("trimap_", ""), kwargs[key]) for key in kwargs.keys() if key.startswith("trimap_")]
        )
        return Trimap(dataset, projection_color, projection_params).plot(), {"display": "none"}, {"display": "none"}, {"display": "block"}, {"display": "none"}

    elif projection_type == "PCA":
        projection_params = dict(
            [(key.replace("pca_", ""), kwargs[key]) for key in kwargs.keys() if key.startswith("pca_")]
        )
        return PCA(dataset, projection_color, projection_params).plot(), {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "block"}

    assert False, "Invalid projection type"

projection_layout = go.Layout(
    paper_bgcolor='white',
    plot_bgcolor='white',
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
)

if __name__ == '__main__':
    app.run_server(debug=True, host="127.0.0.1")