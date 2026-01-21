from typing import Any
from abc import ABC
from markupsafe import Markup

from ...core.templates import catalog

class InputField(ABC):
    def __init__(
        self, name: str, label: str, component: str, type: str,
        id: str | None = None, default: Any = None, 
        required: bool = True,
        pydantic_type: Any = str
    ):
        self.name = name
        self.type = type
        self.label = label
        self.component = component
        self.id = id or name
        self.default = default
        self.data = default
        self.errors: list[str] = []
        self.pydantic_type = pydantic_type
        self.required = required
    
    def render(self, container_class="") -> str:
        return Markup(catalog.render(self.component, field=self, container_class=container_class))