from abc import ABC, abstractmethod

import dash

from cellbro.util.Dataset import Dataset


class DashPage(ABC):
    module: str
    title: str
    path: str
    order: int

    def __init__(self, module: str, title: str, path: str, order: int):
        self.module = module
        self.title = title
        self.path = path
        self.order = order

    @abstractmethod
    def plot(self) -> list:
        ...

    @abstractmethod
    def apply(self):
        ...

    @abstractmethod
    def create_layout(self) -> list:
        ...

    @abstractmethod
    def setup_callbacks(self, dash_app):
        ...

    def create(self):
        dash.register_page(
            self.module,
            title=self.title,
            path=self.path,
            order=self.order,
            layout=self.create_layout(),
        )
