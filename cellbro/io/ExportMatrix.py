from abc import ABC, abstractmethod

import pickle
import pandas as pd

import dash
from dash import Dash, html, dcc, Input, Output, State, ctx

import dash_bootstrap_components as dbc

from .ExportData import ExportData
from ..components.CID import CID
from ..io.FileFormat import FileFormat, CSV, TSV, Pickle

class ExportMatrix(ExportData, ABC):
    def __init__(self, cid: CID, dataset, filename_ph, formats: list[FileFormat] = [CSV, TSV, Pickle]):
        ExportData.__init__(self, cid, dataset, filename_ph, formats=formats)

    @property
    @abstractmethod
    def data(self):
        ...

    @property
    @abstractmethod
    def columns(self):
        ...

    @property
    @abstractmethod
    def index(self):
        ...

    @property
    @abstractmethod
    def features(self):
        ... 

    def export(self, format, feature, filename):
        if filename == "":
            filename = self.filename_ph
        # return params
        if format == ".pkl":
            return dcc.send_bytes(
                src=pickle.dumps(self.data[feature]),
                filename=f"{filename}.pkl",
            )
            # return self.df[feature].to_pickle(f"{feature}.pkl")

        if format == ".csv" or format == ".tsv":
            df = pd.DataFrame(self.data[feature], columns=self.columns, index=self.index)
            return dcc.send_data_frame(
                writer=df.to_csv,
                filename=f"{filename}{format}",
                sep="," if format == ".csv" else "\t",
            )

    def setup_callbacks(self, app):
        _id = self.parent_cid.to_str()
        # POPUP
        output = [
            Output(f"{_id}-modal", "is_open"),
            Output(f"{_id}-export", "data"),
        ]
        inputs = dict(
            open=Input(f"{_id}-open", "n_clicks"),
            close=Input(f"{_id}-close", "n_clicks"),
            export=Input(f"{_id}-apply", "n_clicks"),
            is_open=State(f"{_id}-modal", "is_open"),
            feature=State(f"{_id}-feature-select", "value"),
            format=State(component_id=f"{_id}-format-select", component_property="value"),
            filename=State(component_id=f"{_id}-filename-input", component_property="value"),
        )

        @app.dash_app.callback(
            output=output, inputs=inputs,
            prevent_initial_call=True
        )
        def _(open, close, export, is_open, feature, format, filename):
            if open is None:
                return [False, None]

            if ctx.triggered_id == f"{_id}-apply":
                _file = self.export(format, feature, filename)
                return [False, _file]

            return [not is_open, None]

        # FILE TYPE EXTENSION
        output = Output(f"{_id}-filetype-extension", "children")
        inputs = [Input(f"{_id}-format-select", "value")]

        @app.dash_app.callback(output, inputs)
        def _(value):
            return value

    def _params_layout(self):
        _id = self.parent_cid.to_str()
        return [
            html.Div([
                html.Label("Select Feature to Export:", className="param-label"),
                html.Div([
                    dcc.Dropdown(
                        id=f"{_id}-feature-select",
                        options=self.features,
                        value=self.features[0],
                        placeholder="Select",
                        clearable=False,
                    )
                ], className="param-select")
            ], className="param-row-stacked"),

            html.Div([
                html.Label("Select Format:", className="param-label"),
                html.Div([
                    dcc.Dropdown(
                        id=f"{_id}-format-select",
                        options=dict([(f.ext(), f"{f.name()} ({f.desc()}) [{f.ext()}]") for f in self.formats]),
                        value=self.formats[0].ext(), clearable=False
                    ),
                ], className="param-select")
            ], className="param-row-stacked"),

            html.Div([
                html.Label("Filename:", className="param-label"),
                dbc.InputGroup([
                    dbc.Input(
                        id=f"{_id}-filename-input",
                        value="", type="text",
                        placeholder=self.filename_ph,
                    ),
                    dbc.InputGroupText("", id=f"{_id}-filetype-extension")
                ], className="param-select"),
            ], className="param-row-stacked")
        ]

class ExportLayer(ExportMatrix):
    def __init__(self, cid: CID, dataset, filename_ph, formats: list[FileFormat] = [CSV, TSV, Pickle]):
        ExportMatrix.__init__(self, cid, dataset, filename_ph, formats=formats)


    @property
    def data(self):
        return self.dataset.adata.layers

    @property
    def columns(self):
        return self.dataset.adata.var.index

    @property
    def index(self):
        return self.dataset.adata.obs.index

    @property
    def features(self):
        return list(self.dataset.adata.layers.keys())


class ExportProjection(ExportMatrix):
    def __init__(self, cid: CID, dataset, filename_ph, formats: list[FileFormat] = [CSV, TSV, Pickle]):
        ExportMatrix.__init__(self, cid, dataset, filename_ph, formats=formats)

    @property
    def data(self):
        return self.dataset.adata.obsm

    @property
    def columns(self):
        return None

    @property
    def index(self):
        return self.dataset.adata.obs.index

    @property
    def features(self):
        return list(self.dataset.adata.obsm.keys())
