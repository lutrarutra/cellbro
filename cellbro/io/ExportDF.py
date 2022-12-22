import functools
import os
from abc import ABC, abstractmethod

import pickle
import shutil
import pandas as pd

import dash
from dash import Dash, html, dcc, Input, Output, State, ctx
from dash.exceptions import PreventUpdate

import dash_bootstrap_components as dbc

from cellbro.util.DashAction import DashAction
from cellbro.io.ExportData import ExportData
import cellbro.io.FileFormat as ff

class ExportDF(ExportData):
    def __init__(self, dataset, dataset_attr, id, filename_ph, page_id_prefix, formats: list[ff.FileFormat] = [ff.CSV, ff.TSV, ff.Pickle]):
        ExportData.__init__(self, dataset, id, filename_ph, page_id_prefix, formats=formats)

        if type(dataset_attr) == list:
            self._attr = dataset_attr
        elif type(dataset_attr) == str:
            self._attr = dataset_attr.split(".")
        else:
            assert False, "attr must be a list or a string"

    @property
    def df(self):
        return functools.reduce(lambda obj, attr: getattr(obj, attr), self._attr, self.dataset)

    def export(self, format, features, filename):
        if filename == "":
            filename = self.filename_ph

        if len(features) == 0:
            raise PreventUpdate
        # return params
        if format == ".pkl":
            return dcc.send_bytes(
                src=pickle.dumps(self.df[features]),
                filename=f"{filename}.pkl",
            )
            # return self.df[feature].to_pickle(f"{feature}.pkl")

        if format == ".csv" or format == ".tsv":
            return dcc.send_data_frame(
                writer=self.df[features].to_csv,
                filename=f"{filename}{format}",
                sep="," if format == ".csv" else "\t",
            )
            # self.df[feature].to_csv(f"{feature}{format}", sep="," if format == ".csv" else "\t")

    def setup_callbacks(self, app):
        # POPUP
        output = [
            Output(f"{self._id}-modal", "is_open"),
            Output(f"{self._id}-export", "data"),
        ]
        inputs = dict(
            open=Input(f"{self._id}-open", "n_clicks"),
            close=Input(f"{self._id}-close", "n_clicks"),
            export=Input(f"{self._id}-apply", "n_clicks"),
        )
        state = dict(
            is_open=State(f"{self._id}-modal", "is_open"),
            features=State(f"{self._id}-feature-select", "value"),
            format=State(component_id=f"{self._id}-format-select", component_property="value"),
            filename=State(component_id=f"{self._id}-filename-input", component_property="value"),
        )

        @app.dash_app.callback(
            output=output, inputs=inputs, state=state,
            prevent_initial_call=True
        )
        def _(open, close, export, is_open, features, format, filename):
            if open is None:
                return [False, None]

            if ctx.triggered_id == f"{self._id}-apply":
                return [False, self.export(format, features, filename)]

            return [not is_open, None]

        # FILE TYPE EXTENSION
        output = Output(f"{self._id}-filetype-extension", "children")
        inputs = [Input(f"{self._id}-format-select", "value")]

        @app.dash_app.callback(output, inputs)
        def _(value):
            return value

    def _params_layout(self):
        return [
            html.Div([
                html.Label("Select Feature(s) to Export:", className="param-label"),
                html.Div([
                    dcc.Dropdown(
                        id=f"{self._id}-feature-select",
                        options=self.df.columns,
                        value=self.df.columns,
                        placeholder="Select",
                        multi=True,
                        clearable=False,
                    )
                ], className="param-select")
            ], className="param-row-stacked"),

            html.Div([
                html.Label("Select Format:", className="param-label"),
                html.Div([
                    dcc.Dropdown(
                        id=f"{self._id}-format-select",
                        options=dict([(f.ext(), f"{f.name()} ({f.desc()}) [{f.ext()}]") for f in self.formats]),
                        value=self.formats[0].ext(), clearable=False
                    ),
                ], className="param-select")
            ], className="param-row-stacked"),

            html.Div([
                html.Label("Filename:", className="param-label"),
                dbc.InputGroup([
                    dbc.Input(
                        id=f"{self._id}-filename-input",
                        value="", type="text",
                        placeholder=self.filename_ph,
                    ),
                    dbc.InputGroupText("", id=f"{self._id}-filetype-extension")
                ], className="param-select"),
            ], className="param-row-stacked")
        ]
