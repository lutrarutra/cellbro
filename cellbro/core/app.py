import redis
import rq
import flask

import dash
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import ALL, Dash, Input, Output, State, dcc, html

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
        self.server = flask.Flask(__name__)
        self.redis_url = "redis://localhost:6379"
        self.redis_conn = redis.from_url(self.redis_url)
        self.queue = rq.Queue(connection=self.redis_conn)

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
        self.dataset = Dataset("data/vas.h5ad")

        qc_page = qc.QCPage(self.dataset, self, order=1)
        qc_page.create()

        cells_page = cells.CellsPage(self.dataset, self, order=2)
        cells_page.create()

        de_page = de.DEPage(self.dataset, self, order=3)
        de_page.create()

        pca_page = pca.PCAPage(self.dataset, self, order=4)
        pca_page.create()

        scvi_page = scvi.SCVIPage(self.dataset, self, order=5)
        scvi_page.create()

        # create_genes_page(self.dash_app, self.dataset)

        self.dash_app.layout = html.Div(
            [
                dcc.Location(id="url"),
                html.Div(
                    [
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
                                        )
                                        for page in dash.page_registry.values()
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

    def run(self):
        self.dash_app.run_server(debug=True, host="127.0.0.1")


if __name__ == "__main__":
    app = App()
    app.run()
