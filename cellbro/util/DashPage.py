from abc import ABC, abstractmethod

import dash

import cellbro.util.Components as Components

class DashPage(ABC):
    def __init__(self, module: str, title: str, path: str, order: int):
        self.module = module
        self.title = title
        self.path = path
        self.order = order
        self._id = "home" if path == "/" else path[1:]
        self.actions = dict(
            top_sidebar=Components.HideSidebar(id=f"{self._id}-top-sidebar"),
            bot_sidebar=Components.HideSidebar(id=f"{self._id}-bot-sidebar"),
        )

    @abstractmethod
    def create_layout(self) -> list:
        ...

    def setup_callbacks(self, app):
        for key, action in self.actions.items():
            action.setup_callbacks(app)

    def create(self, app):
        dash.register_page(
            self.module,
            title=self.title,
            path=self.path,
            order=self.order,
            layout=self.create_layout(),
        )
        self.setup_callbacks(app)
