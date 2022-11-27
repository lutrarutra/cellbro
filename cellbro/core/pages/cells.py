import dash
from dash import html, dcc, Input, Output, State
import plotly.express as px
import plotly.graph_objects as go

from cellbro.plotting.Projection import Projection, UMAP, TSNE, PCA, Trimap
from cellbro.plotting.Heatmap import Heatmap
from cellbro.plotting.Violin import Violin

def create_page(dash_app, dataset):
    left_sidebar, main_figure = Projection.create_layout(dataset)
    bottom_left_sidebar, bottom_figure = Heatmap.create_layout(dataset)
    violin_layout = Violin.create_layout(dataset)

    layout = [
        html.Div(id="top", className="top", children=[
            left_sidebar, main_figure, violin_layout
            ]),
        html.Div(id="bottom", className="bottom", children=[
            bottom_left_sidebar, bottom_figure
        ])
    ]

    # Projection
    @dash_app.callback(
        output=Projection.get_callback_outputs(),
        inputs=Projection.get_callback_inputs(),
        state=Projection.get_callback_states()
    )
    def _(submit, projection_color, projection_type, **kwargs):
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

    # Heatmap
    @dash_app.callback(
        output=Heatmap.get_callback_outputs(),
        inputs=Heatmap.get_callback_inputs(),
        state=Heatmap.get_callback_states()
    )
    def _(submit, **kwargs):
        return Heatmap(dataset, kwargs).plot()

    # Violin
    @dash_app.callback(
        output=Violin.get_callback_outputs(),
        inputs=Violin.get_callback_inputs()
    )
    def _(feature, groupby):
        return Violin(dataset).plot(groupby=groupby, feature=feature)

    dash.register_page("pages.cells", path="/cells", order=2, layout=layout)



