from typing import Annotated
from .InputField import InputField

from typing import TypeVar

T = TypeVar("T")

class SelectInput(InputField):
    def __init__(
        self, label: str,
        options: list[tuple[T, str]],
        placeholder: str | None = None,
        default: T | None = None,
        name: str | None = None,
        required: bool = True,
    ):
        super().__init__(
            name=name,
            label=label,
            component="search.Select",
            default=default,
            pydantic_type=Annotated[T, str],
            type="select",
            required=required,
        )
        self.placeholder = placeholder
        self.options = options
        