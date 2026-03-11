from typing import Any, Annotated
from .InputField import InputField

from pydantic import EmailStr, StringConstraints

class StringInputField(InputField):
    def __init__(
        self, label: str,
        placeholder: str | None = None,
        max_length: int | None = None,
        min_length: int | None = None,
        default: str | None = None,
        autocomplete: str | None = None,
        pydantic_type: Any = None,
        required: bool = True,
        name: str | None = None,
        type: str = "text",
    ):
        super().__init__(
            name=name,
            label=label,
            component="inputs.String",
            default=default,
            pydantic_type=pydantic_type or Annotated[str, StringConstraints(max_length=max_length, min_length=min_length)],
            type=type,
            required=required,
        )
        self.placeholder = placeholder
        self.autocomplete = autocomplete


class EmailInputField(StringInputField):
    def __init__(
        self, label: str,
        placeholder: str | None = None,
        max_length: int | None = None,
        min_length: int | None = None,
        default: str | None = None,
        name: str | None = None,
        autocomplete: str | None = "email",
        required: bool = True,
    ):
        super().__init__(
            label=label,
            default=default,
            name=name,
            pydantic_type=Annotated[str, EmailStr, StringConstraints(max_length=max_length, min_length=min_length)],
            autocomplete=autocomplete,
            placeholder=placeholder,
            required=required,
            type="email"
        )


class PasswordInputField(StringInputField):
    def __init__(
        self, label: str,
        placeholder: str | None = None,
        name: str | None = None,
        max_length: int | None = None,
        min_length: int | None = None,
        default: str | None = None,
        autocomplete: str | None = "current-password",
        required: bool = True,
    ):
        super().__init__(
            label=label,
            name=name,
            default=default,
            pydantic_type=Annotated[str, StringConstraints(max_length=max_length, min_length=min_length)],
            autocomplete=autocomplete,
            placeholder=placeholder,
            type="password",
            required=required,
        )