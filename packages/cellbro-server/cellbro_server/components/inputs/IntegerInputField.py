from typing import Any
from .InputField import InputField

from pydantic.types import conint

class IntegerInputField(InputField):
    def __init__(
        self, label: str,
        placeholder: str | None = None,
        default: int | None = None,
        min_value: int | None = None,
        max_value: int | None = None,
        autocomplete: str | None = None,
        pydantic_type: Any = None,
        required: bool = True,
        name: str | None = None,
        type: str = "number",
    ):
        super().__init__(
            name=name,
            label=label,
            component="inputs.Integer",
            default=default,
            pydantic_type=pydantic_type or conint(ge=min_value, le=max_value),
            type=type,
            required=required,
        )
        self.placeholder = placeholder
        self.autocomplete = autocomplete