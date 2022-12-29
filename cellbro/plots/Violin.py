import plotly.graph_objects as go
from dash import Input, Output, State, dcc, html

from ..components.DashFigure import DashFigure
from ..util.DashAction import DashAction
from ..components import components
from ..components.DropDown import DropDown

import scout

default_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=5, b=5, l=5, r=5),
)

class PlotViolin(DashAction):
    def plot(self, feature, groupby):
        fig = scout.ply.violin(self.dataset.adata, y=feature, groupby=groupby, layout=default_layout)
        return [fig]

    def setup_callbacks(self, app):
        output = [
            Output(component_id=f"{self.page_id_prefix}-violin-plot", component_property="figure"),
        ]

        inputs = {
            "feature": Input(component_id=f"{self.page_id_prefix}-violin-feature", component_property="value"),
            "groupby": Input(component_id=f"{self.page_id_prefix}-violin-groupby", component_property="value"),
        }

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(feature, groupby):
            return self.plot(feature, groupby)


class Violin(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.actions.update(
            plot_violin=PlotViolin(dataset, self.page_id_prefix, self.loc_class)
        )

        var_names = sorted(list(self.dataset.adata.var_names))
        other = sorted(list(self.dataset.get_numeric()))
        features = other + var_names
        groupbys = self.dataset.get_categoric()

        self.select_feature = DropDown(
            self.page_id_prefix, id=f"{self.page_id_prefix}-violin-feature",
            options=features, default=features[0],
        )
        self.select_groupby = DropDown(
            self.page_id_prefix, id=f"{self.page_id_prefix}-violin-groupby",
            options=groupbys, default=None, clearable=True
        )
        self.actions.update(self.select_feature.actions)
        self.actions.update(self.select_groupby.actions)

    def create_layout(self) -> list:
        type_params = components.FigureHeaderTab(self.page_id_prefix, tab_label="Type", children=[
            # Features
            html.Div([
                html.Label("Feature"),
                self.select_feature.create_layout(),
                self.select_feature.get_stores(),
                # dcc.Dropdown(
                #     features,
                #     value=list(features.keys())[0],
                #     id=f"{self.page_id_prefix}-violin-feature",
                #     clearable=False,
                # )
            ], className="param-row-stacked"),
            # Groupby select
            html.Div([
                html.Label("Group By"),
                self.select_groupby.create_layout(),
                self.select_groupby.get_stores(),
                # dcc.Dropdown(
                #     groupbys,
                #     value=None,
                #     id=f"{self.page_id_prefix}-violin-groupby",
                #     clearable=True,
                # ),
            ], className="param-row-stacked")
        ])

        figure_params = components.FigureHeader(self.page_id_prefix, tabs=[type_params])


        figure = html.Div(
            children=[
                html.Div(
                    children=figure_params.create_layout(), className="fig-header"
                ),
                html.Div([
                    dcc.Loading(type="circle", children=[
                        html.Div(
                            dcc.Graph(
                                id=f"{self.page_id_prefix}-violin-plot", className="secondary-plot"
                            )
                        )
                    ])
                ], className="secondary-body"),
            ],
            className="secondary",
        )

        return figure

    def get_sidebar_params(self) -> list:
        return []
