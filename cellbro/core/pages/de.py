import dash
from dash import html, dcc, Input, Output, State
from dash.exceptions import PreventUpdate

from cellbro.plotting.DE import DE

def create_page(dash_app, dataset):
    top_sidebar, main, secondary = DE.create_layout(dataset)

    layout = [
        html.Div(id="top", className="top", children=[
            top_sidebar, main, secondary
        ]),
        html.Div(id="bottom", className="bottom", children=[
            # bottom_sidebar, bottom_figure
        ])
    ]

    @dash_app.callback(**DE.get_apply_callbacks())
    def _apply(**kwargs):
        return DE(dataset, kwargs).apply()

    @dash_app.callback(**DE.get_plot_callbacks())
    def _plot(**kwargs):
        if kwargs["groupby"] is None: raise PreventUpdate
        return DE(dataset, kwargs).plot()

    dash.register_page("pages.de", title="DE", path="/de", order=4, layout=layout)



