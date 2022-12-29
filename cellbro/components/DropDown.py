from dash import Output, Input, State, ctx, html, dcc
from dash.exceptions import PreventUpdate

from .DashComponent import DashComponent
from ..util.DashAction import DashAction
import time


class HandleActions(DashAction):
    def __init__(self, id, options, default_value, clearable, multi):
        super().__init__(dataset=None, page_id_prefix=None, loc_class="static")
        self.id = id
        self.options = options
        self.default_value = default_value
        self.clearable = clearable
        self.multi = multi
        self.selected_store_id = dict(id="selected-store", component_id=self.id)
        self.input_store_id = dict(id="input-store", component_id=self.id)

    def _update(self, value):
        if value not in self.options:
            if not (value is None and self.clearable):
                if self.default_value not in self.options:
                    if len(self.options) > 0:
                        return self.options[0]
                    else:
                        return None
                else:
                    return self.default_value

        return value

    def setup_callbacks(self, app):
        @app.dash_app.callback(
            output=[
                Output(self.selected_store_id, "data"),
                Output(self.id, "options"),
                Output(self.id, "value"),
            ],
            inputs=[
                Input(self.input_store_id, "data"),
                Input(self.id, "options"),
                Input(self.id, "value"),
            ],
            state=State(self.selected_store_id, "data"),
        )
        def _(input, options, value, selected):
            if selected is None:
                selected = dict(value=value)

            if ctx.triggered_id == self.id:
                selected = dict(value=value)

            elif ctx.triggered_id == self.input_store_id:
                if input is not None:
                    if "options" in input.keys():
                        self.options = input["options"]
                    if "value" in input.keys():
                        selected = dict(value=input["value"])
            
            if isinstance(selected["value"], list):
                temp = []
                for val in selected["value"]:
                    val = self._update(val)
                    if val is not None or self.clearable:
                        temp.append(val)
                selected["value"] = temp
            else:
                selected["value"] = self._update(selected["value"])


            return selected, self.options, selected["value"]


class DropDown(DashComponent):
    def __init__(self, page_id_prefix, id, options, default, clearable=False, multi=False, placeholder=None, style=None):
        super().__init__(page_id_prefix)
        self.id = id
        self.options = options
        self.default = default
        self.clearable = clearable
        self.multi = multi
        self.placeholder = placeholder
        self.style = style
        self.actions.update(
            {f"{self.id}-handle_actions":HandleActions(self.id, self.options, self.default, self.clearable, self.multi)}
        )

    def get_stores(self):
        default_store = dcc.Store(id=dict(id="selected-store", component_id=self.id), storage_type="local")
        input_store = dcc.Store(id=dict(id="input-store", component_id=self.id), storage_type="local")

        return html.Div([default_store, input_store])

    def create_layout(self):
        return dcc.Dropdown(
            options=self.options, value=self.default,
            id=self.id, clearable=self.clearable,
            multi=self.multi, placeholder=self.placeholder,
            style=self.style
        )
