import dash
from dash import html, dcc, Input, Output, State

from cellbro.plotting.PCA import PCA

def create_page(dash_app, dataset):
    top_sidebar, main_figure, secondary_figure, bottom_sidebar, bottom_figure = PCA.create_layout(dataset)

    layout = [
        html.Div(id="top", className="top", children=[
            top_sidebar, main_figure, secondary_figure
            ]),
        html.Div(id="bottom", className="bottom", children=[
            bottom_sidebar, bottom_figure
        ])
    ]

    @dash_app.callback(
        output=PCA.get_callback_outputs(),
        inputs=PCA.get_callback_inputs(),
    )
    def _(color, pc_x, pc_y, hist_n_pcs, hist_type):
        return PCA(dataset, color, pc_x, pc_y, hist_n_pcs, hist_type).plot()

    dash.register_page("pages.pca", title="PCA", path="/pca", order=4, layout=layout)



