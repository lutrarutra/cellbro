from typing import Literal

import dash
from dash import html, dcc, Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from .DashComponent import DashComponent
from ..util.DashAction import DashAction
from .CID import CID, LocClass


class HideSidebar(DashAction):
    def __init__(self, parent_cid: CID, btn_id: str, side: Literal["left", "right"]):
        super().__init__(parent_cid, dataset=None)
        self.btn_id = btn_id
        self.side = side

    def setup_callbacks(self, app):
        @app.dash_app.callback(
            output=Output(self.parent_cid.to_str(), "className"),
            inputs=[
                Input(self.btn_id, "value"),
                State(self.parent_cid.to_str(), "className")
            ],
        )
        def _(value, cls):
            if value != None:
                if value:
                    return cls + f" {self.side}-open"
                else:
                    return cls.replace(f"{self.side}-open", "")

            raise PreventUpdate

class Sidebar(DashComponent):
    def __init__(
        self, page_id, loc_class, title, params_children,
        create_btn=False, btn_text="Apply"
    ):
        super().__init__(page_id, loc_class, "sidebar")
        self.side: Literal["left", "right"] = "right" if self.loc_class == LocClass.secondary else "left"
        self.row: Literal["top", "bot"] = "bot" if self.loc_class == LocClass.bottom else "top"
        self.create_btn = create_btn
        self.apply_btn_id = f"{self.page_id}-{self.loc_class}-sidebar-apply_btn"
        self.title = title
        self.params_children = params_children
        self.btn_text = btn_text
        self.actions = dict(
            hide_sidebar=HideSidebar(
                parent_cid=CID(self.page_id, self.loc_class, self.type),
                btn_id=f"{self.side}-sidebar-btn",
                side=self.side,
            )
        )

    def create_layout(self):
        if self.create_btn:
            btn_container = dbc.Button(
                self.btn_text,
                color="primary",
                className="mr-1",
                id=self.apply_btn_id,
            )
        else:
            btn_container = []

        return html.Div(
            children=[
                html.Div(
                    [
                        html.H3(self.title),
                    ],
                    className="sidebar-header",
                ),
                dcc.Loading(
                    type="circle", parent_className=f"sidebar-container",
                    children=[
                        html.Div(
                            children=self.params_children,
                            className=f"sidebar-parameters",
                        ),
                        html.Div(
                            children=btn_container,
                            className=f"sidebar-footer",
                        ),
                    ],
                ),
            ],
            className=f"{self.side}-sidebar sidebar {self.row}-sidebar", id=self.cid.to_str()
        )