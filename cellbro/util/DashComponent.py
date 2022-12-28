from abc import ABC, abstractmethod


class DashComponent(ABC):
    def __init__(self, page_id_prefix):
        self.actions = {}
        self.page_id_prefix = page_id_prefix

    @abstractmethod
    def create_layout(self):
        ...

    def setup_callbacks(self, app):
        for action in self.actions.values():
            action.setup_callbacks(app)