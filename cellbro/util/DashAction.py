from abc import ABC, abstractmethod

class DashAction(ABC):
    def __init__(self, dataset, page_id_prefix, loc_class=None):
        self.dataset = dataset
        self.page_id_prefix = page_id_prefix
        self.loc_class = loc_class

    @abstractmethod
    def setup_callbacks(self, app):
        ...