import tkinter as tk
from tkinter import filedialog

from dash import Dash, html, dcc, Input, Output, State
import plotly.express as px
import plotly.graph_objects as go

from cellbro.util.Dataset import Dataset
from cellbro.plotting.Projection import UMAP, TSNE, PCA, Trimap

import scanpy as sc

app = Dash(__name__)
app.enable_dev_tools(debug=True, dev_tools_hot_reload=False)

# root = tk.Tk()
# root.withdraw()
# path = filedialog.askopenfilename()
dataset = Dataset("data/vas.h5ad")

app.layout = html.Div([
    html.Div(children=[
            html.Div(children=[
                html.Label("Projection Type"),
                dcc.Dropdown(["UMAP", "Trimap", "t-SNE", "PCA"], value="UMAP", id="projection-type"),
            ], style={"padding": "10px", "flex": "1"}),
            html.Div(children=[
                html.Label("Color"),
                dcc.Dropdown(dataset.adata.obs_keys(), value=dataset.adata.obs_keys()[0], id="projection-color"),
            ], style={"padding": "10px", "flex": "1"})
    ], style={"padding": "10px", "flex": "1", "display": "flex", "flex-direction": "row"}),

    html.Div(children=[
        dcc.Loading(
            id="loading-projection", type="circle",
            children=[html.Div(dcc.Graph(id="projection-plot"))],
        ),
    ], style={"padding": "10px", "flex": "3"}),
    html.Div(id="color-output")
], style={"display": "flex", "flex-direction": "row"})

@app.callback(
    Output(component_id="projection-plot", component_property="figure"),
    Input(component_id="projection-color", component_property="value"),
    Input(component_id="projection-type", component_property="value"),
)
def update_projection(color, projection_type):
    if projection_type == "UMAP":
        return UMAP(dataset, color).plot()

    elif projection_type == "t-SNE":
        return TSNE(dataset, color).plot()
    
    elif projection_type == "PCA":
        return PCA(dataset, color).plot()

    elif projection_type == "Trimap":
        return Trimap(dataset, color).plot()

    assert False, "Invalid projection type"

projection_layout = go.Layout(
    paper_bgcolor='white',
    plot_bgcolor='white',
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
)


if __name__ == '__main__':
    app.run_server(debug=True, host="127.0.0.1")