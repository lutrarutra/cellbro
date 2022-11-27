import dash
from dash import html, dcc, Input, Output, State
import plotly.express as px
import plotly.graph_objects as go

from cellbro.plotting.QC import QC

def create_page(dash_app, dataset):

    bottom_left_sidebar, bottom_figure = QC.create_layout(dataset)

    layout = [
        html.Div(id="top", children=[
            
            ]),
        html.Div(id="bottom", children=[
            bottom_left_sidebar, bottom_figure
        ])
    ]

    @dash_app.callback(
        output=QC.get_callback_outputs(),
        inputs=QC.get_callback_inputs(),
    )
    def _(submit):
        return QC(dataset).plot()

    dash.register_page("pages.qc", title="QC", path="/qc", order=1, layout=layout)



