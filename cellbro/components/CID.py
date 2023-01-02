from enum import Enum
from typing import Union
from dataclasses import dataclass

class LocClass(Enum):
    main = "main"
    secondary = "secondary"
    bottom = "bottom"
    static = "static"           # Static feature
    dummy = "dummy"             # Dummy component, usually hidden

    @classmethod
    def list(cls):
        return [c.value for c in cls]

    def __str__(self) -> str:
        return self.name

@dataclass
class CID:
    _page_id: str
    _loc_class: LocClass
    _type: str

    def __init__(self, page_id: str, loc_class: Union[LocClass, str], _type: str):
        self._page_id = page_id
        if type(loc_class) == str:
            loc_class = LocClass[loc_class]
        self._loc_class = loc_class
        self._type = _type

    def to_dict(self) -> dict:
        return dict(
            page_id=self._page_id,
            loc_class=self._loc_class.name,
            type=self._type,
        )

    def to_str(self) -> str:
        return f"{self._page_id}-{self._loc_class}-{self._type}"

    @property
    def page_id(self) -> str:
        return self._page_id

    @property
    def loc_class(self) -> LocClass:
        return self._loc_class

    @property
    def type(self) -> str:
        return self._type

    