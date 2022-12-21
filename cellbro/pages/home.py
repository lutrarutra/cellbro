import dash
from dash import Dash, html, dcc, Input, Output, State
import dash_bootstrap_components as dbc

import pandas as pd
import numpy as np

from cellbro.util.DashPage import DashPage
from cellbro.util.DashAction import DashAction
import cellbro.util.Components as Components

import cellbro.io as io

class HomePage(DashPage):
    def __init__(self, dataset, order):
        super().__init__("pages.home", "Home", "/", order)
        self.dataset = dataset
        self.import_popups = dict(
        )
        self.export_popups = dict(
            # export_dataset=IO.ExportDataset(
            #     dataset=self.dataset, id="export-dataset",
            #     filename_ph="adata", formats=[ff.H5AD]
            # ),
            export_cell_features=io.ExportDF(
                dataset=self.dataset, dataset_attr="adata.obs", id="export-cell-features",
                filename_ph="cell_features"
            ),
            export_gene_features=io.ExportDF(
                dataset=self.dataset, dataset_attr="adata.var", id="export-gene-features",
                filename_ph="gene_features"
            ),
            export_layer=io.ExportLayer(
                dataset=self.dataset, id="export-layer",
                filename_ph="layer"
            ),
            export_projection=io.ExportProjection(
                dataset=self.dataset, id="export-projection",
                filename_ph="projection"
            )
        )
        self.actions.update(self.import_popups)
        self.actions.update(self.export_popups)

    def create_layout(self) -> list:
        top_sidebar = Components.create_sidebar(
            id="home-top-sidebar", class_name="top-sidebar",
            title="Home", params_children=[],
            btn_id=None, btn_text=""
        )

        bot_sidebar = Components.create_sidebar(
            id="home-bot-sidebar", class_name="bot-sidebar",
            title="Home", params_children=[],
            btn_id=None, btn_text=""
        )

        obs_divs = []
        for feature in self.dataset.adata.obs.columns:
            obs_divs.append(
                _create_div(feature, Type=str(self.dataset.adata.obs.dtypes[feature]))
            )

        var_divs = []
        for feature in self.dataset.adata.var.columns:
            var_divs.append(
                _create_div(feature, Type=str(self.dataset.adata.var.dtypes[feature]))
            )

        main = html.Div([
            html.Div([
                html.Div([
                    html.H2("Dataset"),
                    html.Div([
                        dbc.Button("Import", id="import-dataset-open", color="primary", disabled=True),
                        dbc.Button("Export", id="export-dataset-open", color="primary", disabled=True),
                    ])
                ], className="title"),
                html.Div([
                    html.Div([
                        html.H3("File:", className="param-label"),
                        html.H4(self.dataset.path, className="value")
                    ], className="col"),
                    html.Div([
                        html.H3("Cells:", className="param-label"),
                        html.H4(self.dataset.adata.shape[0]),
                    ], className="col"),
                    html.Div([
                        html.H3("Genes:", className="param-label"),
                        html.H4(self.dataset.adata.shape[1]),
                    ], className="col"),
                ], className="list")
            ], className="floating-box"),

            html.Div([
                html.Div([
                    html.H2("Cell Features"),
                    html.Div([
                        dbc.Button("Import", id="import-cell-features-open", color="primary"),
                        dbc.Button("Export", id="export-cell-features-open", color="primary"),
                    ])
                ], className="title"),
                html.Div(obs_divs, className="list"),
            ], className="floating-box"),

            html.Div([
                html.Div([
                    html.H2("Gene Features"),
                    html.Div([
                        dbc.Button("Import", id="import-gene-features-open", color="primary"),
                        dbc.Button("Export", id="export-gene-features-open", color="primary"),
                    ])
                ], className="title"),
                html.Div(var_divs, className="list"),
            ], className="floating-box"),

        ], className="home-main")

        layer_divs = [
            _create_div("X (log1p counts)")
        ]
        for layer in self.dataset.adata.layers.keys():
            layer_divs.append(
                _create_div(layer)
            )


        emb_divs = []
        for emb in self.dataset.adata.obsm.keys():
            emb_divs.append(
                _create_div(emb, Dimensions=self.dataset.adata.obsm[emb].shape[1])
            )

        genelist_divs = []
        for genelist in self.dataset.get_gene_lists():
            genelist_divs.append(
                _create_div(genelist)
            )

        bottom = html.Div([
            html.Div([
                html.Div([
                    html.H2("Gene Lists"),
                    html.Div([
                        dbc.Button("Import", id="import-gene-list-open", color="primary"),
                        dbc.Button("Export", id="export-gene-list-open", color="primary"),
                    ])
                ], className="title"),
                html.Div(genelist_divs, className="list"),
            ], className="floating-box"),

            html.Div([
                html.Div([
                    html.H2("Layers"),
                    html.Div([
                        dbc.Button("Import", id="import-layer-open", color="primary"),
                        dbc.Button("Export", id="export-layer-open", color="primary"),
                    ])
                ], className="title"),
                html.Div(layer_divs, className="list"),
            ], className="floating-box"),

            html.Div([
                html.Div([
                    html.H2("Projections"),
                    html.Div([
                        dbc.Button("Import", id="import-projection-open", color="primary"),
                        dbc.Button("Export", id="export-projection-open", color="primary"),
                    ])
                ], className="title"),
                html.Div(emb_divs, className="list"),
            ], className="floating-box"),


        ], className="home-bottom")

        layout =  [
            html.Div(
                id="top",
                className="top",
                children=[top_sidebar, main],
            ),
            html.Div(
                id="bottom", className="bottom", children=[bot_sidebar, bottom]
            ),
        ]

        for popup in self.import_popups.values():
            layout.append(popup.get_layout())

        for popup in self.export_popups.values():
            layout.append(popup.get_layout())

        return layout


def _create_div(name, **kwargs):
    divs = [
        html.Div([
            html.H3("Name:", className="param-label"),
            html.H4(name, className="value")
        ], className="param-row"),
    ]

    for key, value in kwargs.items():
        if type is not None:
            divs.append(html.Div([
                html.H3(f"{key}:", className="param-label"),
                html.H4(value, className="value")
            ], className="param-row"))

    return html.Div(divs, className="col")
