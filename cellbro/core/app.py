import dash
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import ALL, Dash, Input, Output, State, dcc, html

import cellbro.pages.de as de
from cellbro.pages.cells import create_page as create_cells_page
from cellbro.pages.genes import create_page as create_genes_page
from cellbro.pages.pca import create_page as create_pca_page
from cellbro.pages.qc import create_page as create_qc_page
from cellbro.util.Dataset import Dataset

# from cellbro.core.pages.home import create_page as create_home_page


class App:
    def __init__(self):
        self.dash_app = Dash(
            __name__,
            use_pages=True,
            pages_folder="../pages",
            assets_folder="../assets",
            external_stylesheets=[dbc.themes.BOOTSTRAP],
        )
        self.dash_app.enable_dev_tools(debug=True, dev_tools_hot_reload=False)
        self.dataset = Dataset("data/vas.h5ad")

        create_cells_page(self.dash_app, self.dataset)
        create_qc_page(self.dash_app, self.dataset)
        create_genes_page(self.dash_app, self.dataset)
        create_pca_page(self.dash_app, self.dataset)
        # create_de_page(self.dash_app, self.dataset)
        de_page = de.DEPage(self.dataset, self.dash_app)
        de_page.create()

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
