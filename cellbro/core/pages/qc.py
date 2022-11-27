import dash
from dash import html, dcc, Input, Output, State
from dash.exceptions import PreventUpdate
import plotly.express as px
import plotly.graph_objects as go

from cellbro.plotting.QC import QC

def create_page(dash_app, dataset):
    top_sidebar, main_figure, secondary_figure, bottom_sidebar, bottom_figure = QC.create_layout(dataset)

    layout = [
        html.Div(id="top", className="top", children=[
            top_sidebar, main_figure, secondary_figure
            ]),
        html.Div(id="bottom", className="bottom", children=[
            bottom_sidebar, bottom_figure
        ])
    ]

    @dash_app.callback(
        output=QC.get_callback_outputs(),
        inputs=QC.get_callback_inputs(),
    )
    def _(**kwargs):
        return QC(dataset, kwargs).plot()

    outputs, inputs, states = QC.get_filtering_callbacks()
    @dash_app.callback(
        output=outputs,
        inputs=inputs,
        state=states,
    )
    def _(submit, **kwargs):
        return QC(dataset, kwargs).filter(submit)

    dash.register_page("pages.qc", title="QC", path="/qc", order=1, layout=layout)



