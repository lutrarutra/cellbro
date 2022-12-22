from abc import ABC, abstractmethod

class DashFigure(ABC):
    def __init__(self, dataset, page_id_prefix):
        self.dataset = dataset
        self.page_id_prefix = page_id_prefix
        self.actions = {}

    def setup_callbacks(self, app):
        print("Setting up callbacks for", self.__class__.__name__)
        for action in self.actions.values():
            action.setup_callbacks(app)

    @abstractmethod
    def create_layout(self) -> list:
        ...

