import plotly.graph_objects as go
import plotly.express as px
from dash import html, dcc, Input, Output, State
import dash_bootstrap_components as dbc

import scanpy as sc

class Heatmap():
    def __init__(self, dataset, params):
        self.dataset = dataset

    def plot(self):
        fig = px.imshow(self.dataset.adata.X)
        return fig
    
        
    @staticmethod
    def get_callback_outputs():
        return [
            Output(component_id="heatmap", component_property="style")
        ]
        
    # Inputs to Projection
    @staticmethod
    def get_callback_inputs():
        return {
            "submit": Input(component_id="heatmap-submit", component_property="n_clicks"),
        }

    @staticmethod
    def get_callback_states():
        states = {}

        return states

    @staticmethod
    def create_layout(dataset):
        sidebar = html.Div([

        ])

        figure = html.Div(children=[
            html.Div([
                dcc.Loading(
                    id="loading-heatmap", type="circle",
                    children=[html.Div(dcc.Graph(id="heatmap-plot"))],
                )
            ], id="bottom-figure")
            ], id="bottom")

        return sidebar, figure