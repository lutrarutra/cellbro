import dash
from dash import Output, Input, State, ctx, html, dcc, ALL
from dash.exceptions import PreventUpdate

from .DashComponent import DashComponent
from ..util.DashAction import DashAction
from .DashStore import DashStore
from .CID import CID

class DropDown(DashComponent):
    def __init__(
        self, cid: CID, options, default, clearable=False, multi=False, placeholder=None,
        style=None, options_callback=None, update_store_id=None, master_cid=None
    ):
        super().__init__(cid.page_id, cid.loc_class, cid.type)
        self.options = options
        self.default = default
        self.clearable = clearable
        self.multi = multi
        self.placeholder = placeholder
        self._update_store_id = ("update_store-" + "".join(self.cid.type.split("-")[1:])) if update_store_id is None else update_store_id
        self.options_callback = options_callback
        self.static = self.options_callback == None
        self.style = style
        self.children.update(
            default_store=DashStore(self.cid.page_id, self.cid.loc_class, type=f"default_store-{self.cid.type}"),
        )
        self.master_cid = master_cid
        if self.master_cid is not None:
            self.static = True

    @property
    def update_store_id(self):
        assert not self.static, "Cannot update a static dropdown"
        return self._update_store_id

    def setup_callbacks(self, app):
        output = [Output(self.cid.to_dict(), "value")]
        
        inputs = dict(
            _=Input(self.children["default_store"].cid.to_dict(), "modified_timestamp"),
            default_store=State(self.children["default_store"].cid.to_dict(), "data"),
            value=State(self.cid.to_dict(), "value"),
        )
        if self.master_cid is not None:
            output.append(Output(self.cid.to_dict(), "options"))
            inputs.update(dict(master=Input(self.master_cid.to_dict(), "value")))

        @app.dash_app.callback(output=output, inputs=inputs)
        def get_default_value(default_store, value, _, master=None):
            if default_store is None:
                return self.default

            if self.master_cid is not None and ctx.triggered_id == self.master_cid.to_dict():
                options = self.options_callback(master)
                if "value" in default_store.keys():
                    value = default_store["value"]
                if value not in options:
                    value = next(iter(options), None)
                return value, options

            if "value" in default_store.keys() and default_store["value"] == value:
                raise PreventUpdate

            if "value" in default_store.keys():
                value = default_store["value"]
                if self.static:
                    options = self.options
                else:
                    options = self.options_callback()

                if not self.multi:
                    if value not in options:
                        if self.default not in options:
                            if len(options) == 0:
                                if self.clearable:
                                    value = None
                                else:
                                    raise PreventUpdate
                            value = options[0]
                        else:
                            value = self.default
                else:
                    value = [v for v in value if v in options]

                assert value is not None
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

        # if self.master_cid is not None:
        #     @app.dash_app.callback(
        #         output=[
        #             Output(self.cid.to_dict(), "options"),
        #             Output(self.cid.to_dict(), "value"),
        #         ], 
        #         inputs=dict(
        #             value=Input(self.master_cid.to_dict(), "value"),
        #         )
        #     )
        #     def update(value):
        #         options = self.options_callback(value)
        #         print(options)
        #         return options, next(iter(options), None)

        if self.static:
            return

        output = Output(self.cid.to_dict(), "options"),
        inputs = [
            Input("url", "pathname"),
            Input(self.update_store_id, "data"),
        ]
        @app.dash_app.callback(output, inputs)
        def update_options(pathname, _):
            options = self.options_callback()
            assert options is not None
            return options

    def get_stores(self):
        return html.Div([
            self.children["default_store"].create_layout()
        ])

    def create_layout(self):
        return dcc.Dropdown(
            options=self.options, value=self.default,
            id=self.cid.to_dict(), clearable=self.clearable,
            multi=self.multi, placeholder=self.placeholder,
            style=self.style
        )
