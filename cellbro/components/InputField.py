from typing import Literal

import dash
from dash import Output, Input, State, ctx, html, dcc, ALL
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from .DashComponent import DashComponent
from .DashStore import DashStore
from .CID import CID


class InputField(DashComponent):
    def __init__(
        self, cid: CID, default, type: Literal["number", "text"], min=None, max=None,
        clearable=False, placeholder=None, style=None
    ):
        super().__init__(cid.page_id, cid.loc_class, cid.type)
        self.default = default
        self.clearable = clearable
        self.min = min
        self.max = max
        self._type = type
        self.placeholder = placeholder

        self.style = style
        self.children.update(
            default_store=DashStore(self.cid.page_id, self.cid.loc_class, type=f"default_store-{self.cid.type}"),
        )

    @property
    def update_store_id(self):
        assert not self.static, "Cannot update a static dropdown"
        return self._update_store_id

    def setup_callbacks(self, app):
        output = Output(self.cid.to_dict(), "value")

        inputs = dict(
            _=Input(self.children["default_store"].cid.to_dict(), "modified_timestamp"),
            default_store=State(self.children["default_store"].cid.to_dict(), "data"),
            value=State(self.cid.to_dict(), "value"),
        )

        @app.dash_app.callback(output=output, inputs=inputs)
        def get_default_value(default_store, value, _):
            if default_store is None:
                return self.default

            if "value" in default_store.keys() and default_store["value"] == value:
                raise PreventUpdate

            if "value" in default_store.keys():
                value = default_store["value"]
                return value

            raise PreventUpdate

        @app.dash_app.callback(
            output=Output(self.children["default_store"].cid.to_dict(), "data"),
            inputs=dict(
                value=Input(self.cid.to_dict(), "value"),
                default_store=State(self.children["default_store"].cid.to_dict(), "data"),
            )
        )
        def update_default_store(value, default_store):
            if default_store is not None and "value" in default_store.keys() and default_store["value"] == value:
                raise PreventUpdate

            return dict(value=value)

    def get_stores(self):
        return html.Div([
            self.children["default_store"].create_layout()
        ])

    def create_layout(self):
        return dbc.Input(
            id=self.cid.to_dict(),
            value=self.default,
            min=self.min,
            max=self.max,
            step=1,
            className="param-input",
            placeholder=self.placeholder,
            style=self.style,
            type=self._type,
        )
