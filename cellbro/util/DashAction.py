from abc import ABC, abstractmethod

# from ..util.Dataset import Dataset
from ..components.CID import CID, LocClass

class DashAction(ABC):
    def __init__(self, parent_cid: CID, dataset):
        self.parent_cid: CID = parent_cid
        self.dataset = dataset

    @property
    def page_id(self) -> str:
        return self.parent_cid.page_id

    @property
    def loc_class(self) -> LocClass:
        return self.parent_cid.loc_class

    @property
    def type(self) -> str:
        return self.parent_cid.type

    @abstractmethod
    def setup_callbacks(self, app):
        ...