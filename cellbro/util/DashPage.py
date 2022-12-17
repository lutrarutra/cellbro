from abc import ABC, abstractmethod

import dash

import cellbro.util.Components as Components

class DashPage(ABC):
    def __init__(self, module: str, title: str, path: str, order: int):
        self.module = module
        self.title = title
        self.path = path
        self.order = order
        self.actions = dict(
            top_sidebar=Components.HideSidebar(id=f"{path[1:]}-top-sidebar"),
            bot_sidebar=Components.HideSidebar(id=f"{path[1:]}-bot-sidebar"),
        )

    @abstractmethod
    def create_layout(self) -> list:
        ...

    def setup_callbacks(self, app):
        for key, action in self.actions.items():
            print(self.path)
            print(key)
            action.setup_callbacks(app)

    def create(self):
        dash.register_page(
            self.module,
            title=self.title,
            path=self.path,
            order=self.order,
            layout=self.create_layout(),
        )
