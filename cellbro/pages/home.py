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
        super().__init__("pages.home", "Home", "home", order, path="/")
        self.dataset = dataset
        self.import_popups = dict(
        )
        self.export_popups = dict(
            # export_dataset=IO.ExportDataset(
            #     dataset=self.dataset, id=f"{self.id}-export-dataset",
            #     filename_ph="adata", formats=[ff.H5AD]
            # ),
            export_cell_features=io.ExportDF(
                dataset=self.dataset, dataset_attr="adata.obs", id=f"{self.id}-export-cell-features",
                filename_ph="cell_features", page_id_prefix=self.id
            ),
            export_gene_features=io.ExportDF(
                dataset=self.dataset, dataset_attr="adata.var", id=f"{self.id}-export-gene-features",
                filename_ph="gene_features", page_id_prefix=self.id
            ),
            export_layer=io.ExportLayer(
                dataset=self.dataset, id=f"{self.id}-export-layer",
                filename_ph="layer", page_id_prefix=self.id
            ),
            export_projection=io.ExportProjection(
                dataset=self.dataset, id=f"{self.id}-export-projection",
                filename_ph="projection", page_id_prefix=self.id
            )
        )
        self.actions.update(self.import_popups)
        self.actions.update(self.export_popups)

    def create_layout(self) -> list:
        self.sidebars["left_sidebar"] = Components.Sidebar(
            page_id_prefix=self.id, row="top", side="left",
            title="Home", params_children=[],
            apply_btn_id=None, btn_text=""
        )

        self.sidebars["bot_sidebar"] = Components.Sidebar(
            page_id_prefix=self.id, row="bot", side="left",
            title="Home", params_children=[],
            apply_btn_id=None, btn_text=""
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
                        dbc.Button("Import", id=f"{self.id}-import-dataset-open", color="primary", disabled=True),
                        dbc.Button("Export", id=f"{self.id}-export-dataset-open", color="primary", disabled=True),
                    ])
                ], className="title"),
                html.Div([
                    html.Div([
                        html.H3("File:", className=f"{self.id}-param-label"),
                        html.H4(self.dataset.path, className="value")
                    ], className="col"),
                    html.Div([
                        html.H3("Cells:", className=f"{self.id}-param-label"),
                        html.H4(self.dataset.adata.shape[0]),
                    ], className="col"),
                    html.Div([
                        html.H3("Genes:", className=f"{self.id}-param-label"),
                        html.H4(self.dataset.adata.shape[1]),
                    ], className="col"),
                ], className="list")
            ], className="floating-box"),

            html.Div([
                html.Div([
                    html.H2("Cell Features"),
                    html.Div([
                        dbc.Button("Import", id=f"{self.id}-import-cell-features-open", color="primary"),
                        dbc.Button("Export", id=f"{self.id}-export-cell-features-open", color="primary"),
                    ])
                ], className="title"),
                html.Div(obs_divs, className="list"),
            ], className="floating-box"),

            html.Div([
                html.Div([
                    html.H2("Gene Features"),
                    html.Div([
                        dbc.Button("Import", id=f"{self.id}-import-gene-features-open", color="primary"),
                        dbc.Button("Export", id=f"{self.id}-export-gene-features-open", color="primary"),
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
                        dbc.Button("Import", id=f"{self.id}-import-gene-list-open", color="primary"),
                        dbc.Button("Export", id=f"{self.id}-export-gene-list-open", color="primary"),
                    ])
                ], className="title"),
                html.Div(genelist_divs, className="list"),
            ], className="floating-box"),

            html.Div([
                html.Div([
                    html.H2("Layers"),
                    html.Div([
                        dbc.Button("Import", id=f"{self.id}-import-layer-open", color="primary"),
                        dbc.Button("Export", id=f"{self.id}-export-layer-open", color="primary"),
                    ])
                ], className="title"),
                html.Div(layer_divs, className="list"),
            ], className="floating-box"),

            html.Div([
                html.Div([
                    html.H2("Projections"),
                    html.Div([
                        dbc.Button("Import", id=f"{self.id}-import-projection-open", color="primary"),
                        dbc.Button("Export", id=f"{self.id}-export-projection-open", color="primary"),
                    ])
                ], className="title"),
                html.Div(emb_divs, className="list"),
            ], className="floating-box"),


        ], className="home-bottom")

        layout =  [
            html.Div(
                className="top",
                children=[self.sidebars["left_sidebar"].create_layout(), main],
            ),
            html.Div(
                className="bottom", children=[self.sidebars["bot_sidebar"].create_layout(), bottom]
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
