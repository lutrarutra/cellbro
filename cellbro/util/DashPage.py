from abc import ABC, abstractmethod

import dash

from cellbro.util.DashAction import DashAction
from cellbro.util.Dataset import Dataset


class DashPage(ABC):
    def __init__(self, module: str, title: str, path: str, order: int):
        self.module = module
        self.title = title
        self.path = path
        self.order = order
        self.actions = {}

    @abstractmethod
    def create_layout(self) -> list:
        ...

    def setup_callbacks(self, dash_app):
        for key, action in self.actions.items():
            action.setup_callbacks(dash_app)

    def create(self):
        dash.register_page(
            self.module,
            title=self.title,
            path=self.path,
            order=self.order,
            layout=self.create_layout(),
        )
