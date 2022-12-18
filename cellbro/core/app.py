import flask

import dash
import dash_bootstrap_components as dbc

from dash import ALL, Dash, Input, Output, State, dcc, html, ctx
from dash.exceptions import PreventUpdate

import cellbro.pages.cells as cells
import cellbro.pages.de as de
import cellbro.pages.qc as qc
import cellbro.pages.pca as pca
import cellbro.pages.scvi as scvi
from cellbro.pages.genes import create_page as create_genes_page
from cellbro.util.Dataset import Dataset

# from cellbro.core.pages.home import create_page as create_home_page

class App:
    def __init__(self):
        # self.conn = redis_conn
        # self.queue = queue
        self.sidebar_open = False
        
        self.server = flask.Flask(__name__)

        self.dash_app = Dash(
            __name__,
            server=self.server,
            serve_locally=True,
            use_pages=True,
            pages_folder="../pages",
            assets_folder="../assets",
            external_stylesheets=[dbc.themes.BOOTSTRAP],
        )
        self.dash_app.enable_dev_tools(debug=True, dev_tools_hot_reload=False)
        self.dataset = Dataset(
            "/home/lutrarutra/Documents/dev/bioinfo/cellbrowser/data/full.h5ad")

        qc_page = qc.QCPage(self.dataset, order=1)
        qc_page.create(self)

        cells_page = cells.CellsPage(self.dataset, order=2)
        cells_page.create(self)

        de_page = de.DEPage(self.dataset, order=3)
        de_page.create(self)

        pca_page = pca.PCAPage(self.dataset, order=4)
        pca_page.create(self)

        scvi_page = scvi.SCVIPage(self.dataset, order=5)
        scvi_page.create(self)

        # create_genes_page(self.dash_app, self.dataset)

        self.dash_app.layout = html.Div(
            [
                dcc.Location(id="url"),
                html.Div(
                    [
                        dcc.Store(id="genelist-store"),
                        html.Div(
                            id="sidebar-btn-container",
                            children=[dbc.Switch(id="sidebar-btn", value=self.sidebar_open)]
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.Div(
                                            dcc.Link(
                                                page["title"],
                                                href=page["relative_path"],
                                                id=f"{page['name']}_nav_link",
                                                className="nav-link",
                                            ),
                                        ) for page in dash.page_registry.values()
                                    ]
                                )
                            ],
                            id="navbar",
                        ),
                        dash.page_container,
                    ],
                    id="page",
                ),
            ]
        )

        # TABS
        @self.dash_app.callback(
            [
                Output(
                    component_id=f"{page['name']}_nav_link",
                    component_property="className",
                )
                for page in dash.page_registry.values()
            ],
            Input("url", "pathname"),
        )
        def _(pathname):
            return [
                "nav-link active" if pathname == page["relative_path"] else "nav-link"
                for page in dash.page_registry.values()
            ]

        # GENE CARD
        @self.dash_app.callback(
            output=[
                Output("gene-list-dropdown", "value"),
                Output("gene-list-dropdown", "options"),
                Output("genelist-store", "data")
            ],
            inputs=[
                Input("gene-list-dropdown", "value"),
                Input("new-gene-list-button", "n_clicks"),
            ],
            state=[
                State("selected-gene", "children"),
                State("new-gene-list-input", "value"),
            ],
        )
        def _(gene_list, create_new_list, selected_gene, new_gene_list_name):
            if ctx.triggered_id == "new-gene-list-button":
                if create_new_list is None:
                    raise PreventUpdate
                if new_gene_list_name is None:
                    raise PreventUpdate
                if new_gene_list_name in self.dataset.get_gene_lists():
                    raise PreventUpdate

                self.dataset.adata.uns["gene_lists"][new_gene_list_name] = [selected_gene]

                return self.dataset.get_gene_lists(selected_gene), self.dataset.get_gene_lists(), {"new": True}

            res = self.dataset.update_gene_lists(selected_gene, gene_list)
            return res, self.dataset.get_gene_lists(), {"new":False}

        # NEW GENE LIST
        @self.dash_app.callback(
            output=Output("heatmap-selected-genelists", "options"),
            inputs=[Input("genelist-store", "data")]
        )
        def _(genelist_store):
            return self.dataset.get_gene_lists()

    def run(self):
        self.dash_app.run_server(debug=True, host="127.0.0.1")
