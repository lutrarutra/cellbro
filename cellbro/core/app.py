import flask

import dash
import dash_bootstrap_components as dbc

from dash import ALL, Dash, Input, Output, State, dcc, html, ctx
import plotly.graph_objects as go
from dash.exceptions import PreventUpdate
from dash_extensions.enrich import DashProxy, MultiplexerTransform

from ..pages.home import HomePage
from ..pages.QCPage import QCPage
from ..pages.CellsPage import CellsPage
from ..pages.DEPage import DEPage
from ..pages.GSEAPage import GSEAPage
from ..pages.PCAPage import PCAPage
from ..pages.SCVIPage import SCVIPage

from ..util.Dataset import Dataset

from ..components.GeneListComponents import CreateGeneListPopup
from ..components.TermComponents import AddGenesFromTermPopUp

# from cellbro.core.pages.home import create_page as create_home_page

class App:
    def __init__(self, path):
        self.dataset = Dataset(path)

        self.server = flask.Flask(__name__)
        self.dash_app = DashProxy(
            __name__,
            server=self.server,
            serve_locally=True,
            use_pages=True,
            pages_folder="../pages",
            assets_folder="../assets",
            external_stylesheets=[dbc.themes.BOOTSTRAP],
            transforms=[MultiplexerTransform()],
        )
        
        home_page = HomePage(self.dataset, order=0)
        home_page.create()

        qc_page = QCPage(self.dataset, order=1)
        qc_page.create()

        cells_page = CellsPage(self.dataset, order=2)
        cells_page.create()

        de_page = DEPage(self.dataset, order=3)
        de_page.create()

        gsea_page = GSEAPage(self.dataset, order=4)
        gsea_page.create()

        pca_page = PCAPage(self.dataset, order=5)
        pca_page.create()

        # scvi_page = SCVIPage(self.dataset, order=6)
        # scvi_page.create()

        genelist_popup = CreateGeneListPopup(self.dataset)
        term_popup = AddGenesFromTermPopUp(self.dataset)

        # Fixes weird Plotly 'value error'-bug with first plot
        _ = go.Figure(layout=dict(template='plotly'))

        print("Setting up layout...")
        self.dash_app.layout = html.Div([
            dcc.Location(id="url"),
            html.Div([
                dcc.Store(id="genelist-store"),
                dcc.Store(id="term-store"),
                dcc.Store(id="de-store"),
                dcc.Store(id="gsea-store"),
                dcc.Store(id="update_store-projection_type"),
                dcc.Store(id="update_store-color"),
                dcc.Store(id="update_store-groupby"),
                dcc.Store(id="update_store-feature"),
                dcc.Store(id="update_store-categoricals"),
                dcc.Store(id="update_store-cluster_cells_by"),
                dcc.Store(id="update_store-genelists"),
                dcc.Store(id="update_store-genes"),
                
                genelist_popup.create_layout(),
                term_popup.create_layout(),

                html.Div(
                    id="left-sidebar-btn-container",
                    children=[dbc.Switch(id="left-sidebar-btn", value=False)]
                ),
                html.Div(
                    id="right-sidebar-btn-container",
                    children=[dbc.Switch(id="right-sidebar-btn", value=False)]
                ),
                html.Div([
                    html.Div([
                        html.Div(
                            dcc.Link(
                                page["title"],
                                href=page["relative_path"],
                                id=f"{page['name']}_nav_link",
                                className="nav-link",
                            ),
                        ) for page in dash.page_registry.values()
                    ])
                ], id="navbar"),
                dash.page_container,
            ], id="page"),
        ])

        print("Setting up callbacks...")
        home_page.setup_callbacks(self)
        qc_page.setup_callbacks(self)
        cells_page.setup_callbacks(self)
        de_page.setup_callbacks(self)
        gsea_page.setup_callbacks(self)
        pca_page.setup_callbacks(self)
        # scvi_page.setup_callbacks(self)

        genelist_popup.setup_callbacks(self)
        term_popup.setup_callbacks(self)
        

        # TABS
        @self.dash_app.callback(
            [Output(f"{page['name']}_nav_link", "className") for page in dash.page_registry.values()],
            Input("url", "pathname"),
        )
        def _(pathname):
            return [
                "nav-link active" if pathname == page["relative_path"] else "nav-link"
                for page in dash.page_registry.values()
            ]

        print("App Ready!")


    def run(self, debug=False):
        # self.dash_app.enable_dev_tools(debug=True, dev_tools_hot_reload=False)
        self.dash_app.run(debug=debug, host="127.0.0.1", dev_tools_hot_reload=False)
        return
