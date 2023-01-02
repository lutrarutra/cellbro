import plotly.graph_objects as go
from dash import Input, Output, State, dcc, html

from ..components.DashPlot import DashPlot
from ..util.DashAction import DashAction
from ..components import components
from ..components.DropDown import DropDown
from ..components.CID import CID

import scout

default_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=5, b=5, l=5, r=5),
)

class PlotViolin(DashAction):
    def __init__(
        self, parent_cid, dataset,
        select_feature_cid, select_groupby_cid,
    ):
        super().__init__(parent_cid, dataset)
        self.select_feature_cid = select_feature_cid
        self.select_groupby_cid = select_groupby_cid


    def plot(self, feature, groupby):
        fig = scout.ply.violin(self.dataset.adata, y=feature, groupby=groupby, layout=default_layout)
        return fig

    def setup_callbacks(self, app):
        output = Output(self.parent_cid.to_dict(), "figure")
        inputs = dict(
            feature=Input(self.select_feature_cid.to_dict(), "value"),
            groupby=Input(self.select_groupby_cid.to_dict(), "value"),
        )

        @app.dash_app.callback(output=output, inputs=inputs)
        def _(feature, groupby):
            return self.plot(feature, groupby)


class Violin(DashPlot):
    def __init__(self, dataset, page_id, loc_class):
        super().__init__(dataset, page_id, loc_class)

        var_names = sorted(list(self.dataset.adata.var_names))
        other = sorted(list(self.dataset.get_numeric()))
        features = other + var_names
        groupbys = self.dataset.get_categoric()

        self.children.update(
            select_feature=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-feature"),
                options=features, default=features[0],
                options_callback=lambda: sorted(self.dataset.get_numeric()) + sorted(list(self.dataset.adata.var_names))
            ),
            select_groupby=DropDown(
                cid=CID(self.page_id, self.loc_class, "select-groupby"),
                options=groupbys, default=None, clearable=True,
                options_callback=lambda: self.dataset.get_categoric()
            )
        )
        self.actions.update(
            plot_violin=PlotViolin(
                self.cid, dataset,
                select_feature_cid=self.children["select_feature"].cid,
                select_groupby_cid=self.children["select_groupby"].cid,
            )
        )


    def create_layout(self) -> list:
        type_params = components.FigureHeaderTab(self.page_id, self.loc_class, tab_label="Type", children=[
            # Features
            html.Div([
                html.Label("Feature"),
                self.children["select_feature"].create_layout(),
                self.children["select_feature"].get_stores(),
            ], className="param-row-stacked"),
            # Groupby select
            html.Div([
                html.Label("Group By"),
                self.children["select_groupby"].create_layout(),
                self.children["select_groupby"].get_stores(),
            ], className="param-row-stacked")
        ])

        figure_params = components.FigureHeader(self.page_id, self.loc_class, tabs=[type_params])

        figure = html.Div(
            children=[
                html.Div(
                    children=figure_params.create_layout(), className="fig-header"
                ),
                html.Div([
                    dcc.Loading(type="circle", children=[
                        html.Div(
                            dcc.Graph(
                                id=self.cid.to_dict(), className=f"{self.loc_class}-plot"
                            )
                        )
                    ])
                ], className=f"{self.loc_class}-body"),
            ],
            className=self.loc_class.name,
        )

        return figure

    def get_sidebar_params(self) -> list:
        return []
