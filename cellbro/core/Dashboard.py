
from dash import html, dcc
import plotly.express as px
import plotly.graph_objects as go

import scanpy as sc

projection_layout = go.Layout(
    paper_bgcolor='white',
    plot_bgcolor='white',
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=False),
)

class Dashboard():
    def __init__(self, app):
        self.app = app

    def _projection_figure(self):
        umap_fig = px.scatter(
            x=self.app.dataset.adata.obsm["X_umap"][:,0], y=self.app.dataset.adata.obsm["X_umap"][:,1],
            color=self.app.dataset.adata.obs["sample"], color_discrete_sequence=sc.pl.palettes.default_20,
            labels={"x": "UMAP 1", "y": "UMAP 2"}
        )
        umap_fig.update_layout(projection_layout)
        return umap_fig


    def create_layout(self):
        layout = html.Div([
            html.Div(children=[
                    html.Label("Color"),
                    dcc.Dropdown(self.app.dataset.adata.obs_keys()),
            ]),

            html.Div(children=[
                dcc.Graph(figure=self._projection_figure())
            ])
        ], style={"display": "flex", "flex-direction": "row"})

        return layout








