import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, ctx
from dash.exceptions import PreventUpdate

from ..util.DashFigure import DashFigure
from ..util.DashAction import DashAction
from ..util import Components
from ..plots import DE
from ..util.GeneListComponents import SelectGene

import scout

class PlotVolcano(DashAction):
    def plot(self, groupby, reference):
        key = f"rank_genes_{groupby}"

        if reference is None:
            reference = list(self.dataset.adata.uns[key].keys())[0]

        fig = scout.ply.marker_volcano(
            self.dataset.adata.uns[key][reference], layout=DE.figure_layout
        )

        return fig

    def setup_callbacks(self, app):
        outputs = [
            Output(component_id=f"{self.page_id_prefix}-{self.loc_class}-plot", component_property="figure"),
            Output(f"{self.page_id_prefix}-{self.loc_class}-groupby", "options"),
            Output(f"{self.page_id_prefix}-{self.loc_class}-groupby", "value"),
            Output(f"{self.page_id_prefix}-{self.loc_class}-reference", "options"),
            Output(f"{self.page_id_prefix}-{self.loc_class}-reference", "value")
        ]
        
        inputs = dict(
            groupby=Input(
                component_id=f"{self.page_id_prefix}-{self.loc_class}-groupby", component_property="value"
            ),
            reference=Input(
                component_id=f"{self.page_id_prefix}-{self.loc_class}-reference", component_property="value"
            ),
            de_store=Input("de-store", "data")
        )

        @app.dash_app.callback(output=outputs, inputs=inputs)
        def _(groupby, reference, de_store):
            groupby_options = self.dataset.get_rank_genes_groups()
            ref_options = []

            if ctx.triggered_id == "de-store":
                print(de_store)
                if de_store["update"]:
                    groupby = de_store["groupby"]


            if groupby is None:
                groupby = next(iter(groupby_options), None)
            
            # if ctx.triggered_id == f"{self.page_id_prefix}-{self.loc_class}-groupby":
            if groupby is not None:
                ref_options = DE.get_reference_options(self.dataset, groupby)
            
            if reference is None:
                reference = next(iter(ref_options), None)

            if reference is None:
                return dash.no_update, groupby_options, groupby, ref_options, reference

            fig = self.plot(groupby, reference)
            return fig, groupby_options, groupby, ref_options, reference


class DEVolcanoFig(DashFigure):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix, loc_class)
        self.actions.update(
            plot_projection=PlotVolcano(self.dataset, self.page_id_prefix, self.loc_class),
            click_action=SelectGene(dataset=self.dataset, page_id_prefix=self.page_id_prefix, loc_class=self.loc_class),
        )

    def create_layout(self):
        groupby_options = self.dataset.get_rank_genes_groups()
        if len(groupby_options) > 0:
            ref_options = DE.get_reference_options(self.dataset)
        else:
            ref_options = []

        select_ref_tab = Components.FigureHeaderTab(self.page_id_prefix, tab_label="Reference", children=[
            html.Div([
                html.Label("Group By"),
                dcc.Dropdown(
                    options=groupby_options, value=next(iter(groupby_options), None),
                    id=f"{self.page_id_prefix}-{self.loc_class}-groupby", clearable=False,
                ),
            ], className="param-row-stacked"),
            # Projection Hue celect
            html.Div([
                html.Label("Reference"),
                dcc.Dropdown(
                    options=ref_options,
                    value=next(iter(ref_options), None),
                    id=f"{self.page_id_prefix}-{self.loc_class}-reference",
                    clearable=False,
                )
            ], className="param-row-stacked")
        ])

        select_gene_tab = Components.FigureHeaderTab(self.page_id_prefix, tab_label="Gene",
            id=f"{self.page_id_prefix}-{self.loc_class}-genecard", children=[
            Components.create_gene_card(None, self.dataset)
        ])

        fig_header = Components.FigureHeader(self.page_id_prefix, tabs=[select_ref_tab, select_gene_tab])

        figure = html.Div(
            children=[
                html.Div(fig_header.create_layout(), className="fig-header"),

                html.Div(
                    [
                        dcc.Loading(
                            type="circle",
                            children=[
                                html.Div(
                                    dcc.Graph(
                                        id=f"{self.page_id_prefix}-{self.loc_class}-plot", className=f"{self.loc_class}-plot"
                                    )
                                )
                            ],
                        )
                    ],
                    className=f"{self.loc_class}-body",
                ),
            ],
            className=f"{self.loc_class}",
        )
        return figure
