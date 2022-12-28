import flask

import dash
import dash_bootstrap_components as dbc

from dash import ALL, Dash, Input, Output, State, dcc, html, ctx
from dash.exceptions import PreventUpdate

import cellbro.pages.home as home
import cellbro.pages.CellsPage as CellsPage
import cellbro.pages.DEPage as DEPage
import cellbro.pages.QCPage as QCPage
import cellbro.pages.PCAPage as PCAPage
import cellbro.pages.SCVIPage as SCVIPage
import cellbro.pages.GSEAPage as GSEAPage
from cellbro.util.Dataset import Dataset
from ..util.GeneListComponents import CreateGeneListPopup
from ..util.TermComponents import AddGenesFromTermPopUp

# from cellbro.core.pages.home import create_page as create_home_page

class App:
    def __init__(self):
        # self.conn = redis_conn
        # self.queue = queue
        self.left_sidebar_open = False
        self.right_sidebar_open = False
        
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
        self.dataset = Dataset("./data/full.h5ad")

        home_page = home.HomePage(self.dataset, order=0)
        home_page.create()

        qc_page = QCPage.QCPage(self.dataset, order=1)
        qc_page.create()

        cells_page = CellsPage.CellsPage(self.dataset, order=2)
        cells_page.create()

        de_page = DEPage.DEPage(self.dataset, order=3)
        de_page.create()

        gsea_page = GSEAPage.GSEAPage(self.dataset, order=4)
        gsea_page.create()

        pca_page = PCAPage.PCAPage(self.dataset, order=5)
        pca_page.create()

        scvi_page = SCVIPage.SCVIPage(self.dataset, order=6)
        scvi_page.create()

        genelist_popup = CreateGeneListPopup(page_id_prefix=None, dataset=self.dataset)
        term_popup = AddGenesFromTermPopUp(page_id_prefix=None, dataset=self.dataset)

        # create_genes_page(self.dash_app, self.dataset)

        self.dash_app.layout = html.Div(
            [
                dcc.Location(id="url"),
                html.Div(
                    [
                        dcc.Store(id="genelist-store"),
                        dcc.Store(id="term-store"),
                        dcc.Store(id="de-store"),
                        dcc.Store(id="gsea-store"),
                        dcc.Store(id="export-store"),
                        dcc.Store(id="import-store"),
                        genelist_popup.create_layout(),
                        term_popup.create_layout(),
                        html.Div(
                            id="left-sidebar-btn-container",
                            children=[dbc.Switch(id="left-sidebar-btn", value=self.left_sidebar_open)]
                        ),
                        html.Div(
                            id="right-sidebar-btn-container",
                            children=[dbc.Switch(
                                id="right-sidebar-btn", value=self.right_sidebar_open,
                            )]
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

        home_page.setup_callbacks(self)
        qc_page.setup_callbacks(self)
        cells_page.setup_callbacks(self)
        de_page.setup_callbacks(self)
        gsea_page.setup_callbacks(self)
        pca_page.setup_callbacks(self)
        scvi_page.setup_callbacks(self)
        genelist_popup.setup_callbacks(self)
        term_popup.setup_callbacks(self)


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

    def run(self):
        self.dash_app.run_server(debug=True, host="127.0.0.1")
