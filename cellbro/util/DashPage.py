from abc import ABC, abstractmethod

import dash

import cellbro.util.Components as Components
import cellbro.plots.DashFigure as DashFigure

class DashPage(ABC):
    def __init__(self, module: str, title: str, id: str, order: int, path=None):
        self.module = module
        self.title = title
        self.order = order
        self.id = id
        
        if path is None:
            self.path = f"/{id}"
        else:
            self.path = path

        self.actions = dict(
            # left_sidebar=Components.HideSidebar(page_id_prefix=self.id, id=f"{self.id}-left-sidebar", btn_id="left-sidebar-btn"),
            # right_sidebar=Components.HideSidebar(page_id_prefix=self.id, id=f"{self.id}-right-sidebar", btn_id="right-sidebar-btn"),
            # bot_sidebar=Components.HideSidebar(page_id_prefix=self.id, id=f"{self.id}-bot-sidebar", btn_id="left-sidebar-btn"),
        )

        self.plots:dict[str, DashFigure.DashFigure] = {}
        self.sidebars:dict[str, Components.Sidebar] = {}

    @abstractmethod
    def create_layout(self) -> list:
        ...

    def setup_callbacks(self, app):
        for key, action in self.actions.items():
            action.setup_callbacks(app)

        for key, plot in self.plots.items():
            plot.setup_callbacks(app)

        for key, sidebar in self.sidebars.items():
            sidebar.setup_callbacks(app)

    def create(self):
        dash.register_page(
            self.module,
            title=self.title,
            path=self.path,
            order=self.order,
            layout=self.create_layout(),
        )
