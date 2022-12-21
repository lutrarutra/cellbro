import functools, os
from abc import ABC, abstractmethod

import pickle
import shutil
import pandas as pd

import dash
from dash import Dash, html, dcc, Input, Output, State, ctx
from dash.exceptions import PreventUpdate

import dash_bootstrap_components as dbc

from cellbro.util.DashAction import DashAction
import cellbro.io.FileFormat as ff


class ExportData(DashAction, ABC):
    def __init__(self, dataset, id, filename_ph, formats: list[ff.FileFormat] = [ff.CSV, ff.TSV, ff.Pickle]):
        DashAction.__init__(self, dataset)

        self._id = id
        self.filename_ph = filename_ph
        self.formats = formats

    @staticmethod
    def setup_callbacks(self, app):
        ...

    @staticmethod
    def _params_layout(self):
        ...

    def get_layout(self):
        return dbc.Modal(id=f"{self._id}-modal", className="popup", is_open=False, children=[
            dcc.Loading(type="circle", children=[
                dbc.ModalHeader(dbc.ModalTitle("Export")),
                dbc.ModalBody(children=self._params_layout()),
                dbc.ModalFooter(
                    html.Div([
                        dbc.Button(
                            "Export", id=f"{self._id}-apply", color="primary"
                        ),
                        dcc.Download(id=f"{self._id}-export"),
                        dbc.Button(
                            "Cancel", id=f"{self._id}-close", color="danger"
                        )
                    ])
                )
            ])
        ])


# class ExportDataset(ExportData):
#     def __init__(self, dataset, id, filename_ph, formats: list[ff.FileFormat] = [ff.H5AD]):
#         ExportData.__init__(self, dataset, id, filename_ph, formats)

#     def export(self, format, filename):
#         if filename == "":
#             filename = self.filename_ph
#         # return params
#         if format == ".h5ad":
#             os.mkdir("tmp")
#             self.dataset.adata.write_h5ad(f"tmp/_temp_.h5ad", compression=None)

#             shutil.make_archive("_temp_", "zip", "tmp")

#             return dcc.send_file(
#                 path="_temp_.zip",
#                 filename=f"{filename}.zip",
#             )


#     def setup_callbacks(self, app):
#         # POPUP
#         output = [
#             Output(f"{self._id}-modal", "is_open"),
#             Output(f"{self._id}-export", "data"),
#         ]
#         inputs = dict(
#             open=Input(f"{self._id}-open", "n_clicks"),
#             close=Input(f"{self._id}-close", "n_clicks"),
#             export=Input(f"{self._id}-apply", "n_clicks"),
#         )
#         state = dict(
#             is_open=State(f"{self._id}-modal", "is_open"),
#             format=State(component_id=f"{self._id}-format-select", component_property="value"),
#             filename=State(component_id=f"{self._id}-filename-input", component_property="value"),
#         )

#         @app.dash_app.callback(
#             output=output, inputs=inputs, state=state,
#             prevent_initial_call=True
#         )
#         def _(open, close, export, is_open, format, filename):
#             if open is None:
#                 return [False, None]

#             if ctx.triggered_id == f"{self._id}-apply":
#                 _file = self.export(format, filename)
#                 return [False, _file]

#             return [not is_open, None]

#         # FILE TYPE EXTENSION
#         output = Output(f"{self._id}-filetype-extension", "children")
#         inputs = [Input(f"{self._id}-format-select", "value")]

#         @app.dash_app.callback(output, inputs)
#         def _(value):
#             return value

#     def _params_layout(self):
#         return [
#             html.Div([
#                 html.Label("Select Format:", className="param-label"),
#                 html.Div([
#                     dcc.Dropdown(
#                         id=f"{self._id}-format-select",
#                         options=dict([(f.ext(), f"{f.name()} ({f.desc()}) [{f.ext()}]") for f in self.formats]),
#                         value=self.formats[0].ext(), clearable=False
#                     ),
#                 ], className="param-select")
#             ], className="param-row-stacked"),

#             html.Div([
#                 html.Label("Filename:", className="param-label"),
#                 dbc.InputGroup([
#                     dbc.Input(
#                         id=f"{self._id}-filename-input",
#                         value="", type="text",
#                         placeholder=self.filename_ph,
#                     ),
#                     dbc.InputGroupText("", id=f"{self._id}-filetype-extension")
#                 ], className="param-select"),
#             ], className="param-row-stacked")
#         ]
