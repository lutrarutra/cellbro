import dash
from dash import html, dcc, Input, Output, State, ctx
from dash.exceptions import PreventUpdate
import plotly.express as px
import plotly.graph_objects as go
import dash_bootstrap_components as dbc

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
    def _plot(**kwargs):
        return QC(dataset, kwargs).plot()

    outputs, inputs, states = QC.get_filtering_callbacks()
    @dash_app.callback(
        output=outputs,
        inputs=inputs,
        state=states,
    )
    def _filter(submit, **kwargs):
        return QC(dataset, kwargs).filter(submit)

    @dash_app.callback(
        [Output("dispersion-info", "children")],
        [Input("dispersion-plot", "clickData")]
    )
    def _click(clickData):
        if clickData is None: raise PreventUpdate
        return QC.on_hover(clickData["points"][0], dataset)

    @dash_app.callback(
        output=[Output("gene-list-dropdown", "value"), Output("gene-list-dropdown", "options")],
        inputs=[Input("gene-list-dropdown", "value"), Input("new-gene-list-button", "n_clicks")],
        state=[State("selected-gene", "children"), State("new-gene-list-input", "value")]
    )
    def _select_gene_list(gene_list, create_new_list, selected_gene, new_gene_list_name):

        if ctx.triggered_id == "new-gene-list-button":
            if create_new_list is None: raise PreventUpdate
            if new_gene_list_name is None: raise PreventUpdate
            if new_gene_list_name in dataset.get_gene_lists(): raise PreventUpdate
            dataset.adata.uns["gene_lists"][new_gene_list_name] = [selected_gene]
            return dataset.get_gene_lists(selected_gene), dataset.get_gene_lists()

        res = dataset.update_gene_lists(selected_gene, gene_list)
        return res, dataset.get_gene_lists()

    # @dash_app.callback(
    #     output=[Output("gene-list-dropdown", "value")],
    #     inputs=[],
    #     state=[, State("selected-gene", "children")]
    # )
    # def _create_gene_list(submit, value, selected_gene):
    #     if submit is None: raise PreventUpdate
    #     if value is None: raise PreventUpdate
    #     if value in dataset.get_gene_lists(): raise PreventUpdate
    #     dataset.adata.uns["gene_lists"][value] = [selected_gene[6:]]
    #     return [dataset.get_gene_lists()]


    dash.register_page("pages.qc", title="QC", path="/qc", order=1, layout=layout)



