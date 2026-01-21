from typing import Any
from abc import ABC
from markupsafe import Markup

from ...core.templates import render_template

class InputField(ABC):
    def __init__(
        self, name: str, label: str, template: str, type: str,
        id: str | None = None, default: Any = None, 
        required: bool = True,
        pydantic_type: Any = str
    ):
        self.name = name
        self.type = type
        self.label = label
        self.template = template
        self.id = id or name
        self.default = default
        self.data = default
        self.errors: list[str] = []
        self.pydantic_type = pydantic_type
        self.required = required
    
    async def render(self, container_class="") -> str:
        return Markup(await render_template(self.template, field=self, container_class=container_class))