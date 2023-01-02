from abc import ABC, abstractmethod
from typing import Literal, Union

from .CID import CID, LocClass
from ..util.DashAction import DashAction


class DashComponent(ABC):
    def __init__(self, page_id: str, loc_class: Union[LocClass, str], type: str):
        self.actions: dict[str, DashAction] = dict()
        self.children: dict[str, DashComponent] = dict()

        if isinstance(loc_class, str):
            assert loc_class in LocClass.list(), f"Invalid loc_class: {loc_class}"
            loc_class = LocClass[loc_class]

        self._cid = CID(page_id, loc_class, type)

    @classmethod
    def from_cid(cls, cid: CID):
        return cls(cid.type, cid.page_id, cid.loc_class)

    @property
    def cid(self) -> CID:
        return self._cid

    @property
    def type(self) -> str:
        return self._cid.type

    @property
    def page_id(self) -> str:
        return self._cid.page_id
    
    @property
    def loc_class(self) -> LocClass:
        return self._cid.loc_class

    @abstractmethod
    def create_layout(self):
        ...

    def setup_callbacks(self, app):
        for _, action in self.actions.items():
            action.setup_callbacks(app)
        
        for _, component in self.children.items():
            component.setup_callbacks(app)
            
