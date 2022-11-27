import dash
from dash import html, dcc, Input, Output, State
import plotly.express as px
import plotly.graph_objects as go


def create_page(dash_app, dataset):
    layout = [
        html.Div(id="top", children=[
            "Genes"
            ]),
        html.Div(id="bottom", children=[
            "Bottom"
        ])
    ]

    dash.register_page("pages.genes", path="/genes", order=3, layout=layout)