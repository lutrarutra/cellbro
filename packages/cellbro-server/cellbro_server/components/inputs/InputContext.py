from abc import ABC
from markupsafe import Markup

from .InputField import InputField
from ...core.templates import catalog

class InputContext(ABC):
    def __init__(self, component: str):
        self.component = component

    @property
    def fields(self) -> list[InputField]:
        fields = []
        for field_name in dir(self):
            field = getattr(self, field_name)
            if isinstance(field, InputField):
                fields.append(field)

        return fields

    def render(self, container_class: str = "") -> str:
        return Markup(catalog.render(self.component, field=self, container_class=container_class))
