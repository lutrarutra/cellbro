from abc import ABC, abstractmethod

import dash

import cellbro.components.components as components
import cellbro.util.DashAction as DashAction

class DashPage(ABC):
    def __init__(self, module: str, title: str, page_id: str, order: int, path=None):
        self.module = module
        self.title = title
        self.order = order
        self.page_id = page_id
        
        if path is None:
            self.path = f"/{page_id}"
        else:
            self.path = path

        self.actions: dict[str, DashAction.DashAction] = {}

        self.components:dict[str, components.DashComponent] = {}

    @abstractmethod
    def create_layout(self) -> list:
        ...

    def setup_callbacks(self, app):
        for key, sidebar in self.components.items():
            sidebar.setup_callbacks(app)
            
        for key, action in self.actions.items():
            action.setup_callbacks(app)


    def create(self):
        dash.register_page(
            self.module,
            title=self.title,
            path=self.path,
            order=self.order,
            layout=self.create_layout(),
        )
