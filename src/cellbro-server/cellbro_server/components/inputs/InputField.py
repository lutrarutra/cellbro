from typing import Any
from abc import ABC
from markupsafe import Markup

from ...core.templates import catalog

class InputField(ABC):
    def __init__(
        self, label: str, type: str,
        id: str | None = None, default: Any = None, 
        required: bool = True,
        component: str | None = None,
        pydantic_type: Any = str,
        name: str | None = None,
    ):
        self.name = name or label.lower().replace(" ", "_")
        self.type = type
        self.label = label
        self.component = component
        self.id = id or self.name
        self.default = default
        self.data = default
        self.errors: list[str] = []
        self.pydantic_type = pydantic_type
        self.required = required
    
    def render(self, container_class="") -> str:
        if not self.component:
            raise ValueError("Component not specified for InputField")
        return Markup(catalog.render(self.component, field=self, container_class=container_class))